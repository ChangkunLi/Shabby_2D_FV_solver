#include "headers.h"

// Testing full state boundary condition

void free_strem_test_1(){
    string filename  = "city0.gri";
    string inputfile = "PARAM.in.full_state_test";

    mesh_t mesh;
    read_gri(filename, mesh);

    // Initialize state
    MatrixXd U_initial = MatrixXd::Ones(3, mesh.nElem);
    U_initial.row(1) *= 0.453;
    U_initial.row(2) *= 0.769;
    cout << "Free-stream full state boundary condition test:" << endl;
    FVsolver(mesh, filename, inputfile, U_initial, "fullstate_test");
}