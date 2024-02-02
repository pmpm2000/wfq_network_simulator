#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <cmath>
#include <random>
#include <iomanip>
#include <ctime>
#include <chrono>

using namespace std;

class Simulator {
public:
    Simulator(string distribution = "CBR", int queue_size = 100, int packet_limit = 1000, double mi = 0.3)
        : packet_limit(packet_limit), distribution(distribution), queue_size(queue_size), mi(mi) {
        A = {0.0, generate_event(1), generate_event(2)};
        A[0] = min(A[1], A[2]) + generate_event(0);
    }

    void start_simulation() {
        cout << "Simulation in progress..." << endl;

        while (packets_served < packet_limit) {
            int event_type = time_algorithm();
            event_algorithm(event_type);
            double last = simulation_percent;
            simulation_percent = round((static_cast<double>(packets_served) / packet_limit) * 100);
            if (simulation_percent != last) {
                cout << simulation_percent << "%" << endl;
            }
        }

        cout << "Simulation completed." << endl;
        calculate_averages();
        print();
    }

private:
    int packet_limit;
    string distribution;
    int queue_size;
    double mi;
    double clock = 0;
    int packets_served = 0;
    int clients_in_queue = 0;
    double last_event_time = 0;
    double queue_busy = 0;
    int status = 0;
    double spacing_time = 0;
    vector<double> total_lambda = {0,0};
    vector<int> number_of_lambda = {0,0};
    vector<double> total_mi = {0.0,0.0};
    vector<double> number_of_mi = {0.0,0.0};
    int number_of_delays = 0;
    int last_index_in_queue = 0;
    double total_delay = 0;

    vector<int> L = {1, 2}; // packet size 

    vector<double> avg_lambda = {0.0,0.0};
    vector<double> avg_serve_time = {0.0,0.0};

    double avg_queue_time = 0;
    double server_busy = 0;
    vector<double> sum_of_bits = {0,0};
    vector<double> inflow = {0.5,0.25}; //speed, e.g. Mb/s
    vector<double> r = {L[0]*inflow[0],L[1]*inflow[1]}; // perfect weights
    double C = r[0]+r[1]; // perfect serving speed


    vector<double> rho = {0.0,0.0};

    int last_type_of_packet = 0;
    double p_type=1;
    double simulation_percent = 0;
    double avg_server_busy = 0;
    double avg_queue_busy = 0;
    double f1_f2 = 0;
    double f2_f1 = 0;
    vector<double> vsi = {0, 0, 0};
    int dropped_paket = 0;
    vector<vector<vector<double>>> arrive_time;
    vector<double> A;

    void print() const {
        cout << "=======================" << endl;
        cout << "Distribution: " << distribution << endl;
        cout << "Packet sizes: [" << L[0] << ", " << L[1] << "]" << endl;
        cout << "r(i): [" << r[0] << ", " << r[1] << "]" << endl;
        cout << "C: " << C << endl;
        cout << "lambda(i): [" << inflow[0] << ", " << inflow[1] << "]" << endl;
        cout << "Clock: " << clock << endl;
        cout << "Status: " << status << endl;
        cout << "Clients in queue: " << clients_in_queue << endl;
        cout << "Last event time: " << last_event_time << endl;
        cout << "Number of delays: " << number_of_delays << endl;
        cout << "Total delay: " << total_delay << endl;
        cout << "Queue busy: " << avg_queue_busy << endl;
        cout << "Average waiting time: " << avg_queue_time << endl;
        cout << "Server busy: " << avg_server_busy << endl;
        cout << "Packets served: " << packets_served << endl;
        cout << "Average lambda 1: " << avg_lambda[0] << endl;
        cout << "Average lambda 2: " << avg_lambda[1] << endl;
        cout << "Rho: 1 " << rho[0] << endl;
        cout << "Rho: 2 " << rho[1] << endl;
        cout << "Average serve time 1: " << avg_serve_time[0] << endl;
        cout << "Average serve time 2: " << avg_serve_time[1] << endl;
        cout << "Dropped packets: " << dropped_paket << endl;
        cout << "Suma bitów 1: " << sum_of_bits[0] << endl;
        cout << "Suma bitów 2: " << sum_of_bits[1] << endl;
        cout << "Weight of inflow 1: " << f1_f2 <<endl;
        cout << "Weight of inflow 2: " << f2_f1 <<endl;
    }

    int time_algorithm() {
        auto minimum = min_element(A.begin(), A.end());
        clock = *minimum;
        return distance(A.begin(), minimum);
    }

    double calculate_vs(int event_type) {
        double L_val = L[event_type - 1];
        return max(spacing_time, vsi[event_type]) + L_val / r[event_type - 1];
    }

    double generate_event(int event_type) {
        if (event_type > 0) {
            double time_between_clients;
            if (distribution == "Poiss") {
                exponential_distribution<double> distribution_exponential((inflow[event_type - 1]));
                time_between_clients = distribution_exponential(random_engine);
            } else if (distribution == "CBR") {
                time_between_clients = 1/inflow[event_type - 1]; //oblczenie odstępu między klientami
            } else {

                return 0.0;
            }
            total_lambda[event_type-1] += time_between_clients;
            number_of_lambda[event_type-1]++;
            return clock + time_between_clients;

        } else if (event_type == 0) {
            double serving_time;
            if(arrive_time.size()>0){
                double p_type = arrive_time[0][1][1];
                serving_time = L[p_type-1]/C; // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                total_mi[p_type-1] += serving_time;
                sum_of_bits[p_type-1]+=L[p_type-1];
                number_of_mi[p_type-1]++;
            }
            else {
                serving_time =L[last_type_of_packet-1]/C;
                total_mi[last_type_of_packet-1] += serving_time;
                sum_of_bits[last_type_of_packet-1]+=L[last_type_of_packet-1];
                number_of_mi[last_type_of_packet-1]++;
            }
            return clock + serving_time;
        }

        return 0.0;
    }

    void pop_from_queue() {
        arrive_time.erase(arrive_time.begin());
    }
    void event_algorithm(int event_type) {
        clients_in_queue = arrive_time.size();
        queue_busy += (clock - last_event_time) * clients_in_queue;
        server_busy += (clock - last_event_time) * status;

        if (event_type > 0) {
            double vs = calculate_vs(event_type);
            vsi[event_type] = vs;

            if (status == 0) {
                last_type_of_packet = event_type;
                number_of_delays++;
                status = 1;
                A[0] = generate_event(0);
            } else if (status == 1) {
                if (arrive_time.size() < queue_size) {
                    vector<double> v_vs = {vs};
                    vector<double> v_clock = {clock, static_cast<double>(event_type)};
                    vector<vector<double>> v_v = {v_vs, v_clock};
                    arrive_time.push_back(v_v);
                    sort(arrive_time.begin(), arrive_time.end());
                } else {
                    dropped_paket++;
                }
            }

            A[event_type] = generate_event(event_type);
        } else if (event_type == 0) {
            if (arrive_time.size() == 0) {
                A[0] = numeric_limits<double>::infinity();
                status = 0;
            } else if (arrive_time.size() > 0) {
                last_type_of_packet = arrive_time[0][1][1];
                spacing_time = arrive_time[0][0][0];
                number_of_delays++;
                double time_arrived = arrive_time[0][1][0];
                total_delay += clock - time_arrived;
                pop_from_queue();
                A[event_type] = generate_event(event_type);
            }

            packets_served++;
        }

        last_event_time = clock;
    }

    void calculate_averages() {

       avg_serve_time[0] = total_mi[0] / number_of_mi[0];
       avg_serve_time[1] = total_mi[1] / number_of_mi[1];
       avg_lambda[0] = total_lambda[0]/number_of_lambda[0];
       avg_lambda[1] = total_lambda[1]/number_of_lambda[1];
       avg_queue_time = total_delay / number_of_delays;
       if (distribution == "Poiss")
        {
         avg_lambda = 1/(total_lambda / number_of_lambda)*2;
        }
        else
        {
        rho[0] = avg_lambda[0] * avg_serve_time[0];
        rho[1] = avg_lambda[1] * avg_serve_time[1];
        }

        avg_queue_busy = queue_busy / clock;
        avg_server_busy = server_busy / clock;
        f1_f2 = sum_of_bits[0]/(sum_of_bits[0]+sum_of_bits[1]);
        f2_f1 = sum_of_bits[1]/(sum_of_bits[0]+sum_of_bits[1]);
    }

private:
   default_random_engine random_engine;

};


int main() {
    Simulator sim;
    sim.start_simulation();

    return 0;
}

