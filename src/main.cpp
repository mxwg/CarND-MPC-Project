#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include <fstream>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "MPC.h"
#include "json.hpp"

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }

double deg2rad(double x) { return x * pi() / 180; }

double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
    auto found_null = s.find("null");
    auto b1 = s.find_first_of("[");
    auto b2 = s.rfind("}]");
    if (found_null != string::npos) {
        return "";
    } else if (b1 != string::npos && b2 != string::npos) {
        return s.substr(b1, b2 - b1 + 2);
    }
    return "";
}

// Evaluate a polynomial.
double polyeval(Eigen::VectorXd coeffs, double x) {
    double result = 0.0;
    for (int i = 0; i < coeffs.size(); i++) {
        result += coeffs[i] * pow(x, i);
    }
    return result;
}

// Fit a polynomial.
// Adapted from
// https://github.com/JuliaMath/Polynomials.jl/blob/master/src/Polynomials.jl#L676-L716
Eigen::VectorXd polyfit(Eigen::VectorXd xvals, Eigen::VectorXd yvals,
                        int order) {
    assert(xvals.size() == yvals.size());
    assert(order >= 1 && order <= xvals.size() - 1);
    Eigen::MatrixXd A(xvals.size(), order + 1);

    for (int i = 0; i < xvals.size(); i++) {
        A(i, 0) = 1.0;
    }

    for (int j = 0; j < xvals.size(); j++) {
        for (int i = 0; i < order; i++) {
            A(j, i + 1) = A(j, i) * xvals(j);
        }
    }

    auto Q = A.householderQr();
    auto result = Q.solve(yvals);
    return result;
}

double norm(double angle) {
    return atan2(sin(angle), cos(angle));
}


int main(int argc, char *argv[]) {
    uWS::Hub h;

    // get error factors from command line
    if (argc < 1 + 5) {
        std::cout << "Not all factors were supplied!" << std::endl;
        exit(1);
    }
    std::vector<double> factors;
    factors.push_back(atof(argv[1]));
    factors.push_back(atof(argv[2]));
    factors.push_back(atof(argv[3]));
    factors.push_back(atof(argv[4]));
    factors.push_back(atof(argv[5]));

    std::ofstream log;
    log.open("error.csv");
    log << "# cte epsi steer_value throttle_value" << std::endl;

    // MPC is initialized here!
    MPC mpc(factors);

    double prevX, prevY;
    size_t steps = 0;

    h.onMessage(
            [&steps,&prevX, &prevY, &log, &mpc](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                                 uWS::OpCode opCode) {
                // "42" at the start of the message means there's a websocket message event.
                // The 4 signifies a websocket message
                // The 2 signifies a websocket event
                if (steps > 500) {  // exit early
                    exit(0);
                }
                steps++;

                string sdata = string(data).substr(0, length);
                cout << sdata << endl;
                if (sdata.size() > 2 && sdata[0] == '4' && sdata[1] == '2') {
                    string s = hasData(sdata);
                    if (s != "") {
                        auto j = json::parse(s);
                        string event = j[0].get<string>();
                        if (event == "telemetry") {
                            double steeringFactor = -1;
                            // j[1] is the data JSON object
                            vector<double> ptsx = j[1]["ptsx"];
                            vector<double> ptsy = j[1]["ptsy"];
                            double px = j[1]["x"];
                            double py = j[1]["y"];
                            double psi = j[1]["psi"];
                            double v = j[1]["speed"];
                            
                            // predict car position with latency
                            int latency_ms = 100;
                            double delta = j[1]["steering_angle"];
                            delta *= steeringFactor;
                            double a = j[1]["throttle"];
                            a *= deg2rad(25);
                            double dt = latency_ms / 1000.0 / 2;
                            double Lf = 2.67;
                            px = (px + v * cos(psi) * dt);
                            py = (py + v * sin(psi) * dt);
                            psi = (psi + v * delta / Lf * dt);
                            v = (v + a * dt);
                            std::cout << j[1]["x"] << "->" << px << " (x)\n";
                            std::cout << j[1]["y"] << "->" << py << " (y)\n";
                            std::cout << j[1]["speed"] << "->" << v << " (v)\n";

                            // transform points to car coordinates
                            Eigen::VectorXd carX(6), carY(6);
                            vector<double> ptsCarX, ptsCarY;
                            double cp = cos(-psi), sp = sin(-psi);

                            for (size_t i = 0; i < ptsx.size(); ++i) {
                                double x = ptsx[i] - px;
                                double y = ptsy[i] - py;
                                carX[i] = x * cp - y * sp;
                                carY[i] = x * sp + y * cp;
                                ptsCarX.push_back(carX[i]);
                                ptsCarY.push_back(carY[i]);
                            }
                            prevX = px;
                            prevY = py;

                            // fit coefficients
                            auto coeffs = polyfit(carX, carY, 3);
//                            std::cout << coeffs << std::endl;

                            // calculate errors
                            double cte = polyeval(coeffs, 0);
                            double epsi = -atan(coeffs[1]);

                            // calculate steering and throttle
                            double steer_value;
                            double throttle_value;

                            Eigen::VectorXd state(6);
                            state << 0, 0, 0, v, cte, epsi;

                            auto solution = mpc.Solve(state, coeffs);

                            steer_value = solution[0] * steeringFactor / deg2rad(25);
                            throttle_value = solution[1];
                            double cost = solution[2] / 1000;

                            std::cout << "steer: " << steer_value << ", throttle: " << throttle_value
                                      << std::endl;

                            json msgJson;
                            // NOTE: Remember to divide by deg2rad(25) before you send the steering value back.
                            // Otherwise the values will be in between [-deg2rad(25), deg2rad(25] instead of [-1, 1].
                            msgJson["steering_angle"] = steer_value;
                            msgJson["throttle"] = throttle_value;

                            // write values to log file
                            log << cte << " " << epsi << " " << steer_value << " " << throttle_value << " " << cost << std::endl;


                            //Display the MPC predicted trajectory
                            vector<double> mpc_x_vals;
                            vector<double> mpc_y_vals;
                            for (size_t i = 3; i < solution.size(); ++i) {
                                if (i % 2 != 0) {
                                    mpc_x_vals.push_back(solution[i]);
                                } else {
                                    mpc_y_vals.push_back(solution[i]);
                                }
                            }

                            //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
                            // the points in the simulator are connected by a Green line

//                            std::cout << "mpc  vals" << std::endl;
//                            for(size_t i = 0; i < mpc_x_vals.size(); ++i)
//                            {
//                                std::cout << mpc_x_vals[i] << " " << mpc_y_vals[i] << std::endl;
//                            }

                            msgJson["mpc_x"] = mpc_x_vals;
                            msgJson["mpc_y"] = mpc_y_vals;

                            //Display the waypoints/reference line
                            vector<double> next_x_vals;
                            vector<double> next_y_vals;

                            for (size_t i = 0; i < 50; ++i) {
                                next_x_vals.push_back(i);
                                next_y_vals.push_back(polyeval(coeffs, i));
                            }

                            //.. add (x,y) points to list here, points are in reference to the vehicle's coordinate system
                            // the points in the simulator are connected by a Yellow line

                            msgJson["next_x"] = next_x_vals;
                            msgJson["next_y"] = next_y_vals;


                            auto msg = "42[\"steer\"," + msgJson.dump() + "]";
//                            std::cout << msg << std::endl;
                            // Latency
                            // The purpose is to mimic real driving conditions where
                            // the car does actuate the commands instantly.
                            //
                            // Feel free to play around with this value but should be to drive
                            // around the track with 100ms latency.
                            //
                            // NOTE: REMEMBER TO SET THIS TO 100 MILLISECONDS BEFORE
                            // SUBMITTING.
                            this_thread::sleep_for(chrono::milliseconds(latency_ms));
                            ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
                        }
                    } else {
                        // Manual driving
                        std::string msg = "42[\"manual\",{}]";
                        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
                    }
                }
            });

    // We don't need this since we're not using HTTP but if it's removed the
    // program
    // doesn't compile :-(
    h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                       size_t, size_t) {
        const std::string s = "<h1>Hello world!</h1>";
        if (req.getUrl().valueLength == 1) {
            res->end(s.data(), s.length());
        } else {
            // i guess this should be done more gracefully?
            res->end(nullptr, 0);
        }
    });

    h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
        std::cout << "Connected!!!" << std::endl;
    });

    h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                           char *message, size_t length) {
        ws.close();
        std::cout << "Disconnected" << std::endl;
    });

    int port = 4567;
    if (h.listen(port)) {
        std::cout << "Listening to port " << port << std::endl;
    } else {
        std::cerr << "Failed to listen to port" << std::endl;
        return -1;
    }
    h.run();
    log.close();
}
