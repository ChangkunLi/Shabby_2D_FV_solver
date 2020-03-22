#include "headers.h"

// Testing wall boundary condition

void free_strem_test_2(){
    string filename  = "city0.gri";
    string inputfile = "PARAM.in.wall_test";

    mesh_t mesh;
    read_gri(filename, mesh);

    // Initialize state
    MatrixXd U_initial = MatrixXd::Ones(3, mesh.nElem);
    U_initial.row(1) *= 0.0;
    U_initial.row(2) *= 0.0;
    cout << "Free-stream wall boundary condition test:" << endl;
    FVsolver(mesh, filename, inputfile, U_initial, "wall_test");
}