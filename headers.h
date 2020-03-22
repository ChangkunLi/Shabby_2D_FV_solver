#include "iostream"
#include "fstream"
#include "sstream"
#include "string"
#include "vector"
#include "cmath"
#include "Eigen/Dense"
#include "Eigen/Sparse"

using namespace Eigen;
using namespace std;

#define g 9.8
#define rho 1000.0

struct boundary_t {
    int nbfgrp;
    MatrixXi nbface;
    MatrixXi nnode;
    vector<string> title;
    vector<MatrixXi> nodes;
};

struct mesh_t {
    int nNode;
    int Dim;
    MatrixXd Node;
    boundary_t B;
    int nElem;
    string QBasis;
    int QOrder;
    MatrixXi Elem;
};

MatrixXd roe(const MatrixXd& UL, const MatrixXd& UR, double nx, double ny, double& smax);
void flux_test();
void read_gri(const string& filename, mesh_t& mesh);
void edgehash(const MatrixXi& E2N, MatrixXi& IE, MatrixXi & BE);
double tri_area(double dX0, double dY0, double dX1, double dY1, double dX2, double dY2);
void process_gri(const string& filename, MatrixXi& I2E, MatrixXi& B2E, MatrixXd& In, MatrixXd& Bn, MatrixXd& Area, MatrixXi& IE, MatrixXi& BE);
void read_param(const string& inputfile, string& btype, double& CFL,
        string& outputR, string& outputF, string& outputU, double& tmin,
        double& tmax, double& dt_out, int& nIteration, double& nTime);
void FVsolver(const mesh_t& mesh, const string& filename, const string& inputfile, MatrixXd U_initial, const string& file_prefix);     
void free_strem_test_1();   
void free_strem_test_2();