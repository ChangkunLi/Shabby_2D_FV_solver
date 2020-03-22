#include "headers.h"

int main(){
    // test implementation of roe flux
    // flux_test();

    // free-stream / free-stream preservation test
    // free_strem_test_1();
    // free_strem_test_2();

    vector<string> Citys;
    Citys.push_back("city0.gri");
    Citys.push_back("city1.gri");
    Citys.push_back("city2.gri");
    for(int icity=0; icity<Citys.size(); icity++){
        string filename  = Citys[icity];
        string inputfile = "PARAM.in";

        mesh_t mesh;
        read_gri(filename, mesh);

        // Initialize state
        MatrixXd U_initial = MatrixXd::Zero(3, mesh.nElem);
        for(int i=0; i<mesh.nElem; i++){
            int n1 = mesh.Elem(i,0);
            int n2 = mesh.Elem(i,1);
            int n3 = mesh.Elem(i,2);
            double x,y;
            x = (mesh.Node(n1-1,0) + mesh.Node(n2-1,0) + mesh.Node(n3-1,0))/3.0;
            y = (mesh.Node(n1-1,1) + mesh.Node(n2-1,1) + mesh.Node(n3-1,1))/3.0;
            U_initial(0,i) = 1.0 + 0.3*exp(-50.0*((x-1.5)*(x-1.5) + (y-0.7)*(y-0.7)));
        }

        cout << filename << " starts" << endl;
        FVsolver(mesh, filename, inputfile, U_initial, "city" + to_string(icity));
        cout << filename << " ends" << endl;
    }
    return 0;
}
