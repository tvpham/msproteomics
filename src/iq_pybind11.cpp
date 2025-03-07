/*
#########################################################################

Author: Thang V. Pham, t.pham@amsterdamumc.nl

All rights reserved.

Citation:

Pham TV, Henneman AA, Jimenez CR. iq: an R package to estimate relative
protein abundances from ion quantification in DIA-MS-based proteomics,
Bioinformatics 2020 Apr 15;36(8):2611-2613.

R software version: 1.9.11

Adapted for Python, 2025

#########################################################################
*/

#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <string>
#include <utility>
#include <exception>
#include <stdexcept>

/** PYTHON begin **/
#include <Eigen/Dense>
#include <csignal>
#include <iostream>
/** PYTHON end **/

#ifdef _OPENMP
#include <omp.h>
#endif

/** PYTHON begin **/
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

// Global variables
volatile sig_atomic_t stop_signal = 0;

/** PYTHON end **/

using namespace std;

//------------------------------ FILTER -------------------------------------


class utils {
public:
    // count the number of unique values in a large INT array
    static int count_unique(int *vec, size_t nrow) {
        int *m = std::max_element(vec, vec + nrow);
        auto tmp = new vector<bool>(*m + 1, false);  // include 0, for a reasonable size m
        for (size_t i = 0; i < nrow; i++) {
            (*tmp)[vec[i]] = true;
        }
        int total = 0;
        for (auto b : (*tmp)) {
            if (b) {
                total++;
            }
        }
        return total;
    }
};


//------------------------------ MaxLFQ -------------------------------------

using Eigen::FullPivHouseholderQR;
using Eigen::HouseholderQR;
using Eigen::MatrixXd;
using Eigen::VectorXd;

class ion_table {

    static int full_connection;
    static FullPivHouseholderQR<MatrixXd> full_qr;

public:
    static void init(int N) {
        MatrixXd AtA = MatrixXd::Zero(N + 1, N + 1);

        for (int j = 0; j < (N - 1); j++) {
            for (int k = j + 1; k < N; k++) {
                AtA(j, k) = -1.0;
                AtA(k, j) = -1.0;
                AtA(j, j) += 1.0;
                AtA(k, k) += 1.0;
            }
        }

        AtA *= 2.0;

        for (int j = 0; j < N; j++) {
            AtA(j, N) = 1.0;
            AtA(N, j) = 1.0;
        }
        AtA(N, N) = 0.0;

        full_qr = FullPivHouseholderQR<MatrixXd>(AtA);

        full_connection = N * (N - 1) / 2;
    }

    unordered_map<int, vector<double>> map;
    int ncol;

    ion_table(int ncol) : ncol(ncol) {
    }

    // add an entry to a map, increase size if necessary
    void add_to_table(int ion, int col, double quant) {
        auto r = map.emplace(ion, vector<double>());
        if (r.second) {
            //n++;
            r.first->second.resize(ncol, NAN);
        }
        r.first->second[col] = quant;
    }

    void spread(int *g, int i, int val) {
        g[i] = val;

        for (const auto &m : map) {
            if (!isnan(m.second[i])) {
                for (int j = 0; j < ncol; j++) {
                    if (g[j] < 0 && !isnan(m.second[j])) {
                        spread(g, j, val);
                    }
                }
            }
        }
    }

    inline double get_median(double *median_buffer, int median_size) {
        sort(median_buffer, median_buffer + median_size);
        int mid = median_size / 2;
        return median_size % 2 == 0 ? (median_buffer[mid] + median_buffer[mid - 1]) / 2.0 : median_buffer[mid];
    }

    void maxLFQ(double *buffer, int *g) {
        if (map.size() < 1) {
            fill_n(buffer, ncol, NAN);
            fill_n(g, ncol, 0);
            return;
        }

        if (map.size() == 1) {
            for (int i = 0; i < ncol; i++) {
                buffer[i] = map.begin()->second[i];
            }
            fill_n(g, ncol, 0);
            return;
        }

        double *median_buffer = new double[map.size()];
        int median_size = 0;

        fill_n(g, ncol, -1);

        int val = 0;
        for (int i = 0; i < ncol; i++) {
            if (g[i] < 0) {
                spread(g, i, val++);
            }
        }

        fill_n(buffer, ncol, NAN);

        for (int i = 0; i < val; i++) {
            vector<int> ind;
            for (int j = 0; j < ncol; j++) {
                if (g[j] == i) {
                    ind.push_back(j);
                }
            }

            if (ind.size() == 1) {
                median_size = 0;

                for (const auto &m : map) {
                    if (!isnan(m.second[ind[0]])) {
                        median_buffer[median_size++] = m.second[ind[0]];
                    }
                }
                if (median_size == 0) {
                    buffer[ind[0]] = NAN;
                } else {
                    buffer[ind[0]] = get_median(median_buffer, median_size);
                }
            } else {

                int N = ind.size();

                MatrixXd AtA = MatrixXd::Zero(N + 1, N + 1);
                VectorXd Atb = VectorXd::Zero(N + 1);

                int n_connection = 0;

                for (int j = 0; j < (N - 1); j++) {
                    for (int k = j + 1; k < N; k++) {
                        median_size = 0;
                        for (const auto &m : map) {
                            if (!isnan(m.second[ind[j]]) && !isnan(m.second[ind[k]])) {
                                median_buffer[median_size++] = m.second[ind[k]] - m.second[ind[j]];
                            }
                        }
                        if (median_size > 0) {
                            n_connection++;

                            double r_i_j = get_median(median_buffer, median_size);
                            AtA(j, k) = -1.0;
                            AtA(k, j) = -1.0;

                            AtA(j, j) += 1.0;
                            AtA(k, k) += 1.0;

                            Atb(j) -= r_i_j;
                            Atb(k) += r_i_j;
                        }
                    }
                }

                AtA *= 2.0;
                Atb *= 2.0;

                for (int j = 0; j < N; j++) {
                    AtA(j, N) = 1.0;
                    AtA(N, j) = 1.0;
                }
                AtA(N, N) = 0.0;

                // mean data
                double sum = 0.0;
                int count = 0;
                for (const auto &m : map) {
                    for (const auto &s : ind) {
                        if (!isnan(m.second[s])) {
                            sum += m.second[s];
                            count++;
                        }
                    }
                }
                Atb(N) = sum * (double)N / (double)count;

                if (n_connection == full_connection) {
                    VectorXd x = full_qr.solve(Atb);
                    for (int j = 0; j < N; j++) {
                        buffer[ind[j]] = x(j);
                    }
                } else {
                    FullPivHouseholderQR<MatrixXd> qr(AtA);
                    VectorXd x = qr.solve(Atb);
                    for (int j = 0; j < N; j++) {
                        buffer[ind[j]] = x(j);
                    }
                }
            }
        }

        delete[] median_buffer;
    }
};

int ion_table::full_connection;
FullPivHouseholderQR<MatrixXd> ion_table::full_qr;


/** PYTHON begin **/

// Signal handler to catch interrupts
void signal_handler(int signum) {
    if (signum == SIGINT) {
        stop_signal = 1;
    }
}

// Check for user interruption
bool check_interrupt() {
    return stop_signal != 0;
}

py::dict iq_MaxLFQ(py::array_t<int> &vp,
                   py::array_t<int> &vi,
                   py::array_t<int> &vs,
                   py::array_t<double> &vq) {
    
    py::gil_scoped_acquire acquire;

    int stop_sig = 0;

    int* proteins = (int*)vp.request().ptr;
    int* ions = (int*)vi.request().ptr;
    int* samples = (int*)vs.request().ptr;
    double* quants = (double*)vq.request().ptr;

    size_t nrow = vp.size();

    int n_proteins = utils::count_unique(proteins, nrow);
    int n_samples = utils::count_unique(samples, nrow);

    vector<double> buffer(n_proteins * n_samples);

    auto group_annotation = new vector<string>(n_proteins, "");
    int* row_names = new int[n_proteins];
    int* col_names = new int[n_samples];

    //py::print("nrow = ", nrow, "; # proteins = ", n_proteins,"; # samples = ", n_samples);

    auto protein_index = new vector<vector<int>>(*(std::max_element(proteins, proteins + nrow)));  // allowing for missing proteins

    for (auto i : *protein_index) {
        i.clear();
    }

    vector<int> sample_set;  // sample index, in case a subset of samples is quantified

    for (size_t i = 0; i < nrow; i++) {
        (*protein_index)[proteins[i] - 1].push_back(i);

        if (samples[i] >= (int)sample_set.size()) {
            sample_set.resize(samples[i] + 1, -1);
            sample_set[samples[i]] = 0;
        }
        else {
            sample_set[samples[i]] = 0;
        }
    }

    int cc = 0;
    for (int i = 0; i < (int)sample_set.size(); i++) {
        if (sample_set[i] > -1) {
            col_names[cc] = i;
            sample_set[i] = cc++;
        }
    }

    ion_table::init(n_samples);  // QR for full matrix

    int nr = 0;
    vector<int> map_back((*protein_index).size(), -1);
    for (size_t i = 0; i < (*protein_index).size(); i++) {
        if (!(*protein_index)[i].empty()) {
            map_back[i] = nr;
            row_names[nr++] = i + 1;
        }
    }

    std::signal(SIGINT, signal_handler);


    size_t thres_display = 0;

    #ifdef _OPENMP

        omp_set_dynamic(0);

        int no_threads = omp_get_num_procs() - 1;

        if (no_threads < 1) {
            no_threads = 1;
        }

        omp_set_num_threads(no_threads);

        //py::print("Using", no_threads, "threads...\n");

    //#else

        //py::print("Using a single thread...\n");

    #endif


    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < (*protein_index).size(); i++) {
        if (stop_signal) {
            continue;
        }
        
        int thread_id = 0;

        #ifdef _OPENMP
            thread_id = omp_get_thread_num();
        #endif

        if (!(*protein_index)[i].empty()) {

            ion_table *tab = new ion_table(n_samples);

            double *w = new double[n_samples];
            int *gr = new int[n_samples];

            for (auto j : (*protein_index)[i]) {
                tab->add_to_table(ions[j], sample_set[samples[j]], quants[j]);
            }

            tab->maxLFQ(w, gr);

            //--- annotation
            bool single_compoment = true;
            int component = -1;

            for (int j = 0; j < n_samples; j++) {
                if (!isnan(w[j])) {
                    if (component == -1) {
                        component = gr[j];
                    } else {
                        if (component != gr[j]) {
                            single_compoment = false;
                            break;
                        }
                    }
                }
            }

            if (component == -1) {
                (*group_annotation)[map_back[i]] = "NA";
            } else {
                if (!single_compoment) {  // all NOT in the same component, otherwise the default empty string is good
                    for (int j = 0; j < n_samples; j++) {
                        if (j > 0) {
                            (*group_annotation)[map_back[i]].push_back(';');
                        }
                        if (isnan(w[j])) {
                            (*group_annotation)[map_back[i]].append("NA");
                        } else {
                            (*group_annotation)[map_back[i]].append(to_string(gr[j] + 1));
                        }
                    }
                }
            }

            //--- filling buffer
            for (int j = 0; j < n_samples; j++) {  
                //size_t k = j * n_proteins + map_back[i];  // R by column
                size_t k = map_back[i]*n_samples + j;       // Python by row
                if (isnan(w[j])) {
                    buffer[k] = NAN;
                } else {
                    buffer[k] = w[j];
                }
            }

            delete[] w;
            delete[] gr;
            delete tab;
        }

        if (thread_id == 0) {
            if (check_interrupt()) {  // user interrupted ...
                stop_signal = 1;
                #pragma omp flush(stop_signal)
            }
        }
    }
    
    if (stop_signal) {
        py::print("Canceled.\n");

        delete protein_index;
        delete group_annotation;
        delete[] col_names;
        delete[] row_names;

       return (py::none());
    }

    py::array mat = py::array(buffer.size(), buffer.data());
    py::array _r_names = py::array(n_proteins, row_names);
    py::array _c_names = py::array(n_samples, col_names);

    py::list ann(n_proteins);
    for (int i = 0; i < n_proteins; i++) {
        ann[i] = group_annotation->at(i);
    }

    py::dict vec;
    vec["estimate"] = mat;
    vec["annotation"] = ann;
    vec["_r_names"] = _r_names;
    vec["_c_names"] = _c_names;

    delete protein_index;
    delete group_annotation;
    delete[] col_names;
    delete[] row_names;
    
    return (vec);
}

extern "C"
{
    #include "main.c"
}

void siteloc(py::list inlist) {

    int argc = (int)inlist.size();

    void *ptr = malloc(argc * sizeof(char*));
    
    if(ptr == NULL) {
        py::print("Cannot allocate memory.\n");
        return;
    }

    char** argv = (char**)ptr;

    for (int i = 0; i < argc; ++i) {
        argv[i] = (char*)PyUnicode_AsUTF8(inlist[i].ptr());
    }

    int ret = main(argc, argv);

    free(ptr);
}



PYBIND11_MODULE(_msproteomics, m) {        
    m.def("iq_maxLFQ", &iq_MaxLFQ, "MaxLFQ in C++");
    m.def("siteloc", &siteloc, "siteloc in C++");
}

/** PYTHON END **/