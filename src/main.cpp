#include <cstdlib>
#include "LBM.h"
#include "IOobject.h"

using namespace std;

int main() {
    // Dictionary for parameters
    map<string, string> input;

    // Read the parameters from input file
    ReadData r("./input.dat");
    r.read(input);

    // Parameters
    int nx = stoi(input["nx"]);
    int ny = stoi(input["ny"]);
    int timeStep = stoi(input["timeStep"]);
    int Re = stoi(input["Re"]);
    int solver = stoi(input["solver"]);
    double utop = stod(input["utop"]);
    double omegaP = stod(input["omegaP"]);
    double omega3 = stod(input["omega3"]);
    double omega4 = stod(input["omega4"]);
    double nu = utop * nx / Re;
    double omega = 1.0 / (3.0 * nu + 0.5);

    // Initialize the LBM solver
    LBM cavity(ny, nx, utop, omega);

    // Run the LBM solver
    switch (solver) {
        case 1:
            cout << "Solver used: BGK Raw Moments" << endl << endl;
            cavity.runBGKRaw(timeStep);
            break;
        case 2:
            cout << "Solver used: MRT Central Moments" << endl << endl;
            cavity.runMRTCentral(timeStep, omegaP, omega3, omega4);
            break;
        case 3:
            cout << "Solver used: BGK Central Moments" << endl << endl;
            cavity.runBGKCentral(timeStep);
            break;
        default:
            cout << "Choose a correct solver:" << endl
                << "1 - BGK Raw Moments" << endl
                << "2 - MRT Central Moments" << endl
                << "3 - BGK Central Moments" << endl;
            break;
    }

    // Write the velocities
    cavity.writeUx();
    cavity.writeUy();

    return 0;
}
