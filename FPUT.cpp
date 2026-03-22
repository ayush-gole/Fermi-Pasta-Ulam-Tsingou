#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

vector<double> acc(vector<double>& q, double alpha){
    int size = q.size();
    vector<double> A_(size, 0.0);
    for (int i = 1; i < size - 1; i++){
        double prev_q = q[i-1];
        double now_q = q[i];
        double next_q = q[i + 1];

        A_[i] = (next_q - 2 * now_q + prev_q) + alpha * ((next_q - now_q) * (next_q - now_q) - (now_q - prev_q) * (now_q - prev_q));
    }
    return A_;
}

vector<double> get_E(int N, vector<double>& q, vector<double>& P, vector<double>& w, vector<vector<double>>& phi, vector<double>& Q_t, vector<double>& P_t){
    
    vector<double> E_tot;
    Q_t.clear();
    P_t.clear();

    for (int k = 1; k < N + 1; k++){
        
        double Q_ = 0, p_ = 0;

        for ( int j = 1; j < N+1; j++){
            Q_ += phi[j-1][k-1] * q[j];
            p_ += phi[j-1][k-1] * P[j];
        }
        double E_K = 0.5 * ( p_ * p_ + w[k] * w[k]  * Q_ * Q_);
        E_tot.push_back(E_K);
        Q_t.push_back(Q_);
        P_t.push_back(p_);

    }

    return E_tot;

}

int main(){
    // Declaring The veriables 
    int N, m, maxT;
    float alpha, dt, A;
    // declaring the constant PI
    const double PI = 3.14159265358979323846;
    // Declaring the arrays for future
    vector<double>time, q, vel, q_prev, TE_arr, C_time, ET1, ET2, ET3, w;
    vector<vector<double>>phi, norm_E;
    
    N = 32; // Number of points 
    m = 1; // Mass of each point is 1 unit
    alpha = 0.25; // A constant for verlet integration

    dt = 0.05; // declaring the Time step 
    maxT = 50000; // Maximum units of Time to be counted
    A = 0.1; // A Constant for Future

    for ( double t = 0; t < maxT; t += dt){
        time.push_back(t);
    }

    q.resize(N + 2, 0);
    vel.resize(N + 2, 0);
    q_prev.resize(N + 2, 0);
    cout<<q.size()<<endl;

    for(int i = 1; i <= N; i++){
        q[i] = A * sin((PI * i) / (N + 1));
        q_prev[i] = q[i];
    }

    phi.resize(N, vector<double>(N, 0.0));
    for(int j = 1;  j < N; j++){
        for(int k = 1; k < N; k++){
            phi[j-1][k-1] = sqrt(2.0 / (N + 1)) * sin((PI * k * j) / (N + 1));
        }
    }
    cout<<phi.size() << endl;

    w.resize(N + 1);
    w[0] = 0.0;
    for(int k = 1; k <= N; k++){
        w[k] = 2.0 * sin((PI * k) / ( 2.0 * (N + 1)));
    }
    cout<<w.size() << endl;

    vector<double> acc_prev = acc(q, alpha);

    // Get initial Q and P values (at t=0)
    vector<double> Q0, P0;
    get_E(N, q, vel, w, phi, Q0, P0);
    // cout<<Q0<<","<<P0;

    vector<double> R_time;
    vector<vector<double>>RK_time;

    for( int i = 0; i < time.size(); i++){
        double t = i * dt;

        for( int j = 1; j < q.size()-1; j++){
            q[j] += vel[j] * dt + 0.5 * acc_prev[j] * dt * dt;
        }

        vector<double> new_acc = acc(q, alpha);
        for( int j = 1; j < vel.size()-1; j++){
            vel[j] += 0.5 * dt * ( acc_prev[j] + new_acc[j]);
        }

        acc_prev = new_acc;

        if (i % 100 == 0){
            C_time.push_back(t);
            vector<double> Q1, P1;
            vector<double> EK = get_E(N, q, vel, w, phi, Q1, P1);

            ET1.push_back(EK[0]);
            ET2.push_back(EK[1]);
            ET3.push_back(EK[2]);

            double tot_EK = 0;
            for(double energy : EK){
                tot_EK += energy;
            }
            vector<double> norm_E_(EK.size(), 0.0);
            for(int i = 0; i<EK.size(); i++){
                norm_E_[i] = EK[i] / tot_EK;
            }

            double nume = 0;
            double deno = 0;
            vector<double> RKT;
            for(int k = 0; k < Q1.size(); k++){
                nume += Q1[k] * Q0[k] + P1[k] * P0[k];
                deno += Q0[k] * Q0[k] + P0[k] * P0[k];
                RKT.push_back(nume/deno);
            }
            RK_time.push_back(RKT);

            double R_t = nume / deno;
            R_time.push_back(R_t);

            norm_E.push_back(norm_E_);

            TE_arr.push_back(tot_EK);
        }
    }

    string name = "FPUT_" + to_string(A) + ".csv";
    ofstream file(name);
    file << "C_time,R_time,ET1,ET2,ET3"; // header row

    for(int mode = 1; mode <= N; mode++){
        file << ",norm_E" << mode;
    }
    for(int mode = 1; mode <= N; mode++){
        file << ",RK_t" << mode;
    }

    file << "\n";

    for (int i = 0; i < C_time.size(); i++) {
        file << C_time[i] <<"," << R_time[i] <<"," << ET1[i] <<"," << ET2[i] <<"," << ET3[i];
        for(int mode = 0; mode < N; mode++){
            file << "," << norm_E[i][mode];
        }
        for(int mode = 0; mode < N; mode++){
            file << "," << RK_time[i][mode];
        }
        file << "\n";
    }
    file.close();

    // ofstream heatmap("heatmap.csv");
    // heatmap << "time_index,mode,normalized_energy\n";
    // for(int i = 0; i < C_time.size(); i++){
    //     for(int mode = 0; mode < N; mode++){
    //         heatmap << i << "," << (mode + 1) << "," << norm_E[i][mode] << "\n";
    //     }
    // }
    // heatmap.close();

    cout<<"Data saved to "<< name << endl;
    return 0;
}