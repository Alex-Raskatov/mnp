#include <ios>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <string>

const int T = 1000;
const int K_X = 100, K_Y = 100;
const int V_X = 5, V_Y = 5;
const int N_VX = 2*V_X + 1, N_VY = 2*V_Y + 1;

const double dvx = 0.1, dvy = 0.1;
const int step = 10;
const int r = 50;
const int V0_X = -3, V0_Y = 0;

const int wall_x_left = 34, wall_x_right = 50;
const int wall_bottom_y = 35, wall_top_y = 65;

int index(int k_x, int k_y, int i_vx, int i_vy) {
    return k_x + k_y*K_X + (i_vx + V_X)*K_X*K_Y + (i_vy + V_Y)*K_X*K_Y*N_VX;
}

double normal_f(int i_vx, int i_vy) {
    return std::exp(-i_vx*i_vx / 2.0) * std::exp(-i_vy*i_vy/2.0) / 10.0;
}

void set_initial_values(double* data) {
    for (int i_vx = -V_X; i_vx < V_X + 1; i_vx++) {
        for (int i_vy = -V_Y; i_vy < V_Y + 1; i_vy++) {
            for (int k_x = 0; k_x < K_X; k_x++) {
                for (int k_y = 0; k_y < K_Y; k_y++) {
                    if (k_x < wall_x_left) {
                        data[index(k_x, k_y, i_vx, i_vy)] = normal_f(i_vx, i_vy);
                    } else {
                        data[index(k_x, k_y, i_vx, i_vy)] = 0;
                    }
                }
            }
        
        }
    }
}

double get_denom_x() {
    double denom = 0.0;
    for (int i_vy = -V_Y; i_vy < V_Y + 1; i_vy++) {
        for (int i_vx = 1; i_vx < V_X + 1; i_vx++) {
            double v = dvx * i_vx;
            denom += v * std::exp(-v*v/2);
        }
    }
    return denom;
}

double get_denom_y() {
    double denom = 0.0;
    for (int i_vx = -V_X; i_vx < V_X + 1; i_vx++) {
        for (int i_vy = 1; i_vy < V_Y + 1; i_vy++) {
            double v = dvy * i_vy;
            denom += v * std::exp(-v*v/2);
        }
    }
    return denom;
}

double get_nom_x(double *current, int k_x, int k_y, int sign) {
    double nom = 0;
    for (int i_vy = -V_Y; i_vy < V_Y + 1; i_vy++) {
        for (int i_vx = -V_X; i_vx < V_X + 1; i_vx++) {
            double v = i_vx * dvx;
            if (sign * v > 0) {
                nom += v * current[index(k_x, k_y, i_vx, i_vy)];
            }
        }
    }
    return sign*nom;
}

double get_nom_y(double *current, int k_x, int k_y, int sign) {
    double nom = 0;
    for (int i_vy = -V_Y; i_vy < V_Y + 1; i_vy++) {
        for (int i_vx = -V_X; i_vx < V_X + 1; i_vx++) {
            double v = i_vy * dvy;
            if (sign * v > 0) {
                nom += v * current[index(k_x, k_y, i_vx, i_vy)];
            }
        }
    }
    //std::cout << "flag3" << '\n';
    return sign*nom;
}


void make_iteration_x (double* current, double* next) {
    for (int i_vx = -V_X; i_vx < V_X + 1; i_vx++) {
        for (int i_vy = -V_Y; i_vy < V_Y + 1; i_vy++) {
            double v = i_vx*dvx;
            if (v > 0) {
                for (int k_x = 0; k_x < K_X; k_x++) {
                    for (int k_y = 0; k_y < K_Y; k_y++) {
                        if ((k_x == wall_x_left or k_x == wall_x_right) and (k_y <= wall_bottom_y or k_y >= wall_top_y) or k_x == 0) {
                            //next[index(k_x, k_y, i_vx, i_vy)] = next[index(k_x, k_y, -i_vx, i_vy)];
                            next[index(k_x, k_y, i_vx, i_vy)] = get_nom_x(current, k_x, k_y, -1)/get_denom_x()*std::exp(-v*v/2);
                        } else {
                            next[index(k_x, k_y, i_vx, i_vy)] = current[index(k_x, k_y, i_vx, i_vy)] - v*(current[index(k_x, k_y, i_vx, i_vy)] - current[index(k_x - 1, k_y, i_vx, i_vy)]);
                        }
                    }
                }
            } else {
                for (int k_x = 0; k_x < K_X; k_x++) {
                    for (int k_y = 0; k_y < K_Y; k_y++) {
                        if ((k_x == wall_x_right or k_x == wall_x_left) and (k_y <= wall_bottom_y or k_y >= wall_top_y) /*or k_x == K_X - 1*/) {
                            //next[index(k_x, k_y, i_vx, i_vy)] = next[index(k_x, k_y, -i_vx, i_vy)];
                            next[index(k_x, k_y, i_vx, i_vy)] = get_nom_x(current, k_x, k_y, 1)/get_denom_x()*std::exp(-v*v/2);
                        } else if (k_x == K_X-1) {
                            next[index(k_x, k_y, i_vx, i_vy)] = 0;
                        } else {
                            next[index(k_x, k_y, i_vx, i_vy)] = current[index(k_x, k_y, i_vx, i_vy)] - v*(current[index(k_x + 1, k_y, i_vx, i_vy)] - current[index(k_x, k_y, i_vx, i_vy)]);
                        }
                    }
                }
            }
        }
    }
}

void make_iteration_y (double* current, double* next) {
    for (int i_vx = -V_X; i_vx < V_X + 1; i_vx++) {
        for (int i_vy = -V_Y; i_vy < V_Y + 1; i_vy++) {
            double v = i_vy*dvy;
            if (v > 0) {
                for (int k_x = 0; k_x < K_X; k_x++) {
                    for (int k_y = 0; k_y < K_Y; k_y++) {
                        if ((k_y == 0 and k_x < wall_x_left) or (k_x >= wall_x_left and k_x <= wall_x_right and (k_y == wall_bottom_y or k_y == wall_top_y))) {
                            //next[index(k_x, k_y, i_vx, i_vy)] = next[index(k_x, k_y, i_vx, -i_vy)];
                            //std::cout << "flag" << '\n';
                            next[index(k_x, k_y, i_vx, i_vy)] = get_nom_y(current, k_x, k_y, -1)/get_denom_y()*std::exp(-v*v/2);
                            //std::cout << "flag" << '\n';
                        } else if (k_y == 0 and k_x >= wall_x_left) {
                            next[index(k_x, k_y, i_vx, i_vy)] = 0;
                        } else {
                            next[index(k_x, k_y, i_vx, i_vy)] = current[index(k_x, k_y, i_vx, i_vy)] - v*(current[index(k_x, k_y, i_vx, i_vy)] - current[index(k_x, k_y - 1, i_vx, i_vy)]);
                        }
                    }
                }
            } else {
                for (int k_x = 0; k_x < K_X; k_x++) {
                    for (int k_y = 0; k_y < K_Y; k_y++) {
                        if ((k_y == K_Y - 1 and k_x < wall_x_left) or (k_x >= wall_x_left and k_x <= wall_x_right and (k_y == wall_top_y or k_y == wall_bottom_y))) {
                            //next[index(k_x, k_y, i_vx, i_vy)] = next[index(k_x, k_y, i_vx, -i_vy)];
                            //std::cout << "flag2" << '\n';
                            next[index(k_x, k_y, i_vx, i_vy)] = get_nom_y(current, k_x, k_y, 1)/get_denom_y()*std::exp(-v*v/2);
                            //std::cout << "flag2" << '\n';
                        } else if (k_y == K_Y - 1 and k_x >= wall_x_left) {
                            next[index(k_x, k_y, i_vx, i_vy)] = 0;
                        } else {
                            next[index(k_x, k_y, i_vx, i_vy)] = current[index(k_x, k_y, i_vx, i_vy)] - v*(current[index(k_x, k_y + 1, i_vx, i_vy)] - current[index(k_x, k_y, i_vx, i_vy)]);
                        }
                    }
                }
            }
        }
    }
}

void save_to_file(double* data, int i) {

    char filename[30];
    sprintf(filename, "data/out_%03d.dat", i);
    std::ofstream file(filename);
    for (int k_x = 0; k_x < K_X; k_x++) {
        for (int k_y = 0; k_y < K_Y; k_y++) {
            double concetration = 0.0;
            for (int i_vx = -V_X; i_vx < V_X + 1; i_vx++) {
                for (int i_vy = -V_Y; i_vy < V_Y + 1; i_vy++) {
                    concetration += data[index(k_x, k_y, i_vx, i_vy)];
                }
            }
            file << k_x << " " << k_y << " " << concetration << std::endl;
        }
    }
    file.close();
}

void save_to_file_number(double* data, int i) {

    std::ofstream file("data_other/number.txt", std::ios::app);
    double number = 0.0;
    for (int k_x = 0; k_x < K_X; k_x++) {
        for (int k_y = 0; k_y < K_Y; k_y++) {
            for (int i_vx = -V_X; i_vx < V_X + 1; i_vx++) {
                for (int i_vy = -V_Y; i_vy < V_Y + 1; i_vy++) {
                    number += data[index(k_x, k_y, i_vx, i_vy)];
                }
            }
        }
    }
    file << number << std::endl;
    file.close();
}

void save_to_file_kin(double* data, int i) {

    char filename[30];
    sprintf(filename, "data_kin/out_%03d.dat", i);
    std::ofstream file(filename);
    for (int k_x = 0; k_x < K_X; k_x++) {
        for (int k_y = 0; k_y < K_Y; k_y++) {
            double energy = 0.0;
            for (int i_vx = -V_X; i_vx < V_X + 1; i_vx++) {
                for (int i_vy = -V_Y; i_vy < V_Y + 1; i_vy++) {
                    energy += data[index(k_x, k_y, i_vx, i_vy)]*((i_vx*dvx)*(i_vx*dvx) + (i_vy*dvy)*(i_vy*dvy))/2;
                }
            }
            file << k_x << " " << k_y << " " << energy << std::endl;
        }
    }
    file.close();
}



int main() {

    double* current = new double[K_X * K_Y * N_VX * N_VY];
    double* next = new double[K_X * K_Y * N_VX * N_VY];

    set_initial_values(current);

    for (int i = 0; i < T; i++) {

        if (i % step == 0) {
            save_to_file(current, i);
            save_to_file_number(current, i);
            save_to_file_kin(current, i);
            std::cout << "Шаг " << i << '\n';
        }
        if (i % 2 == 0) {
            make_iteration_x(current, next);

        } else {
            make_iteration_y(current, next);
        }
        std::swap(current, next);
    }

    delete [] current;
    delete [] next;

    return 0;
}