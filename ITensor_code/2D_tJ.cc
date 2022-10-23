#include "itensor/all.h"
#include <cmath>
#include <vector>
#include <complex>
#include <iostream>
#include <fstream>
#include <sstream>

#define M_PI 3.14159265358979323846 // pi

using namespace itensor;
using namespace std;

int main(int argc, char *argv[])
{
    // H = -t1 Cdag_i C_i+1 + -t2 Cdag_i C_i+2 + j1 (S_i dot S_i+1 - 1/4 n_i n_+1) + j2 (S_i dot S_i+2 - 1/4 n_i n_+2)
    //   + jk S_i dot (S_j cross S_k)
    // Sweeps
    // Parse the input file
    //
    if (argc < 2)
    {
        printfln("Usage: %s SweepTable2Dj1j2", argv[0]);
        return 0;
    }
    auto input = InputGroup(argv[1], "input");

    // Read in individual parameters from the input file
    auto Nx = input.getInt("Nx");
    auto Ny = input.getInt("Ny");
    // Define variables
    double j1 = 1.0;                                   // nn exchange
    auto t1 = input.getReal("t1");                     // nearest neighbor hopping
    auto j2_ini = input.getReal("j2");                     // nnn exchange
    auto t2 = input.getReal("t2", sqrt(j2_ini / j1) * t1); // next nearest neighbor hopping
    auto j2_scale = input.getReal("j2scale", 1);       // 0 for no j2
    auto j2 = j2_ini * j2_scale;
    cout << "j1 = " << j1 << "j2 = " << j2 << " t1 = " << t1 << " t2 = " << t2 << endl;
    auto jk = input.getReal("jk");              // chiral interaction
    auto doping = input.getInt("DopingNumber"); // electron density = (N - doping)/N
    if (doping % 2 == 1)
    {
        Error("doping must be even number");
    }
    auto hz = input.getReal("hz");               // edge magnetic pining field
    auto jkedge = input.getReal("jkedge");       // edge Chiral pining field
    auto torus = input.getYesNo("torus", false); // torus
    int N = Nx * Ny;
    auto nsweeps = input.getInt("nsweeps");
    auto quiet = input.getYesNo("quiet", true);
    auto yperiodic = input.getYesNo("periodic", true);
    auto lattice_type = input.getString("lattice_type", "triangularXC"); // triangularXC
    auto ErrorGoal = input.getReal("ErrorGoal");
    auto ReadPsiSites = input.getYesNo("ReadPsiSitesFromFile");
    auto DoDMRG = input.getYesNo("DoDMRG", false);
    auto DoMidDMRGSweeps = input.getYesNo("DoMidDMRGSweeps", false);
    auto SavePsiSites = input.getYesNo("SavePsiSites", true);
    auto DoChiralityTriangle = input.getYesNo("DoChiralityTriangle", false);
    auto DoSpinStructure = input.getYesNo("DoSpinStructure", false);
    auto DoCorrelation = input.getYesNo("DoCorrelation", false);
    auto TypeofCorrelation = input.getString("TypeofCorrelation", "S"); // S, Sz, Sxy
    auto SCcorr = input.getYesNo("SCcorrelation", false);               // different types
    auto DoVEE = input.getYesNo("DoVEE", false);
    auto DoMidEnergyPerSite = input.getYesNo("DoMidEnergyPerSite", false);
    auto DoLocalSpinValue = input.getYesNo("DoLocalSpinValue", false);
    auto CmtIn = input.getString("CommentsInPsiSites", "testin");
    auto CmtOut = input.getString("CommentsOutPsiSitesCorr", "testout");
    bool Writememorytodisk = true;
    int whentowrite = 6000;

    auto table = InputGroup(input, "sweeps");
    // Create the sweeps class & print
    //
    auto sweeps = Sweeps(nsweeps, table);
    println(sweeps);
    // set variables to string type
    //
    std::ostringstream strNxs;
    strNxs << Nx;
    string strNx = strNxs.str();
    std::ostringstream strNys;
    strNys << Ny;
    string strNy = strNys.str();
    std::ostringstream strj2s;
    strj2s << j2_ini;
    std::string strj2 = strj2s.str();
    std::ostringstream strj2_scales;
    strj2_scales << j2_scale;
    std::string strj2_scale = strj2_scales.str();
    std::ostringstream strdopings;
    strdopings << doping;
    std::string strdoping = strdopings.str();
    std::ostringstream strjks;
    strjks << jk;
    std::string strjk = strjks.str();
    std::ostringstream strt1s;
    strt1s << t1;
    std::string strt1 = strt1s.str();

    // Define space
    //
    // tJ sites
    tJ sites;
    if (ReadPsiSites)
    {
        cout << "Reading sites 0.5 from file" << endl;
        readFromFile("sites_tJ_" + lattice_type + "_Nx_" + strNx + "_Ny_" + strNy + "_doping_" + strdoping + "_t1_" + strt1 + "_j2_" + strj2 + "_j2scale_" + strj2_scale + "_jk_" + strjk + CmtIn, sites);
    }
    else
    {
        sites = tJ(N, {"ConserveQNs=", true});
    }

    Args args;
    args.add("Quiet", quiet);
    args.add("YPeriodic", yperiodic);
    args.add("maxdim", 500);
    args.add("Cutoff", 1E-8);
    println(args);

    Args args1;
    args1.add("Quiet", quiet);
    if (Writememorytodisk)
    {
        args1.add("WriteDim", whentowrite);
    }

    // Initialize the site degrees of freedom.
    //
    cout << "total site number: " << length(sites) << endl;
    cout << "MPO" << endl;
    auto ampo = AutoMPO(sites);
    auto ampo_mid = AutoMPO(sites);
    int mid_from = Nx / 4 + 1; // Nx/4+1
    int mid_to = Nx * 3 / 4;   // Nx*3/4
    if (Nx % 4 == 0)
    {
        mid_from = Nx / 4 + 1;
        mid_to = Nx * 3 / 4;
    }
    else if (Nx % 2 == 0)
    {
        mid_from = (Nx + 2) / 4 + 1;
        mid_to = (Nx + 2) * 3 / 4 - 1;
    }
    else
    {
        Error("Nx must be even number");
    }
    // specific value
    // mid_from =
    // mid_to =

    if (lattice_type == "triangularXC")
    {
        //
        //       /\  /\  /\  /
        //      /__\/__\/__\/
        //     /\  /\  /\  /
        //    /__\/__\/__\/
        //
        auto lattice = triangularLattice(Nx, Ny, args);
        // nn interaction
        for (int n = 1; n <= N; n++)
        {
            int x = (n - 1) / Ny + 1;
            int y = (n - 1) % Ny + 1;
            vector<vector<int>> nn(2);
            if (x == 1)
            {
                if (y == Ny)
                {
                    if (yperiodic)
                    {
                        nn[0] = {n, n};
                        nn[1] = {n - Ny + 1, n + Ny};
                    }
                    else
                    {
                        nn[0] = {n};
                        nn[1] = {n + Ny};
                    }
                }
                else
                {
                    nn[0] = {n, n};
                    nn[1] = {n + 1, n + Ny};
                }
            }
            else if (x == Nx)
            {
                if (y == Ny)
                {
                    if (yperiodic)
                    {
                        nn[0] = {n, n};
                        nn[1] = {n - Ny + 1, n - 2 * Ny + 1};
                    }
                    else
                    {
                        nn[0] = {};
                        nn[1] = {};
                    }
                }
                else
                {
                    nn[0] = {n, n};
                    nn[1] = {n + 1, n - Ny + 1};
                }
            }
            else
            {
                if (y == Ny)
                {
                    if (yperiodic)
                    {
                        nn[0] = {n, n, n};
                        nn[1] = {n - Ny + 1, n + Ny, n - 2 * Ny + 1};
                    }
                    else
                    {
                        nn[0] = {n};
                        nn[1] = {n + Ny};
                    }
                }
                else
                {
                    nn[0] = {n, n, n};
                    nn[1] = {n + 1, n + Ny, n - Ny + 1};
                }
            }
            for (int i = 0; i < nn[0].size(); i++)
            {
                ampo += (0.5 * j1), "S+", nn[0][i], "S-", nn[1][i];
                ampo += (0.5 * j1), "S-", nn[0][i], "S+", nn[1][i];
                ampo += (-0.5 * j1), "Nup", nn[0][i], "Ndn", nn[1][i];
                ampo += (-0.5 * j1), "Ndn", nn[0][i], "Nup", nn[1][i];
                ampo += -t1, "Cdagup", nn[0][i], "Cup", nn[1][i];
                ampo += -t1, "Cdagup", nn[1][i], "Cup", nn[0][i];
                ampo += -t1, "Cdagdn", nn[0][i], "Cdn", nn[1][i];
                ampo += -t1, "Cdagdn", nn[1][i], "Cdn", nn[0][i];

                if ((DoMidEnergyPerSite) && ((nn[0][i] - 1) / Ny + 1 > mid_from - 1) && ((nn[0][i] - 1) / Ny + 1 <= mid_to) && ((nn[1][i] - 1) / Ny + 1 > mid_from) && ((nn[1][i] - 1) / Ny + 1 <= mid_to))
                {
                    ampo_mid += (0.5 * j1), "S+", nn[0][i], "S-", nn[1][i];
                    ampo_mid += (0.5 * j1), "S-", nn[0][i], "S+", nn[1][i];
                    ampo_mid += (-0.5 * j1), "Nup", nn[0][i], "Ndn", nn[1][i];
                    ampo_mid += (-0.5 * j1), "Ndn", nn[0][i], "Nup", nn[1][i];
                    ampo_mid += -t1, "Cdagup", nn[0][i], "Cup", nn[1][i];
                    ampo_mid += -t1, "Cdagup", nn[1][i], "Cup", nn[0][i];
                    ampo_mid += -t1, "Cdagdn", nn[0][i], "Cdn", nn[1][i];
                    ampo_mid += -t1, "Cdagdn", nn[1][i], "Cdn", nn[0][i];
                }
            }
        }
        // nnn interaction for triangular lattice
        if ((j2 != 0) || (t2 != 0))
        {
            for (int n = 1; n <= N; n++)
            {
                int x = (n - 1) / Ny + 1;
                int y = (n - 1) % Ny + 1;
                vector<vector<int>> nn(2);
                if (x == 1)
                {
                    if (y == Ny)
                    {
                        if (yperiodic)
                        {
                            nn[0] = {n};
                            nn[1] = {n + 1};
                        }
                        else
                        {
                            nn[0] = {};
                            nn[1] = {};
                        }
                    }
                    else
                    {
                        nn[0] = {n};
                        nn[1] = {n + Ny + 1};
                    }
                }
                else if (x == 2)
                {
                    if (y == Ny)
                    {
                        if (yperiodic)
                        {
                            nn[0] = {n, n};
                            nn[1] = {n + 1, n - 2 * Ny + 2};
                        }
                        else
                        {
                            nn[0] = {};
                            nn[1] = {};
                        }
                    }
                    else if (y == (Ny - 1))
                    {
                        if (yperiodic)
                        {
                            nn[0] = {n, n};
                            nn[1] = {n + Ny + 1, n - 2 * Ny + 2};
                        }
                        else
                        {
                            nn[0] = {n};
                            nn[1] = {n + Ny + 1};
                        }
                    }
                    else
                    {
                        nn[0] = {n, n};
                        nn[1] = {n + Ny + 1, n - Ny + 2};
                    }
                }
                else if (x == Nx)
                {
                    if (y == Ny)
                    {
                        if (yperiodic)
                        {
                            nn[0] = {n, n};
                            nn[1] = {n - 2 * Ny + 2, n - 3 * Ny + 1};
                        }
                        else
                        {
                            nn[0] = {};
                            nn[1] = {};
                        }
                    }
                    else if (y == (Ny - 1))
                    {
                        if (yperiodic)
                        {
                            nn[0] = {n, n};
                            nn[1] = {n - 2 * Ny + 2, n - 2 * Ny + 1};
                        }
                        else
                        {
                            nn[0] = {n};
                            nn[1] = {n - 2 * Ny + 1};
                        }
                    }
                    else
                    {
                        nn[0] = {n, n};
                        nn[1] = {n - 2 * Ny + 1, n - Ny + 2};
                    }
                }
                else
                {
                    if (y == Ny)
                    {
                        if (yperiodic)
                        {
                            nn[0] = {n, n, n};
                            nn[1] = {n - 2 * Ny + 2, n - 3 * Ny + 1, n + 1};
                        }
                        else
                        {
                            nn[0] = {};
                            nn[1] = {};
                        }
                    }
                    else if (y == (Ny - 1))
                    {
                        if (yperiodic)
                        {
                            nn[0] = {n, n, n};
                            nn[1] = {n - 2 * Ny + 2, n - 2 * Ny + 1, n + Ny + 1};
                        }
                        else
                        {
                            nn[0] = {n, n};
                            nn[1] = {n - 2 * Ny + 1, n + Ny + 1};
                        }
                    }
                    else
                    {
                        nn[0] = {n, n, n};
                        nn[1] = {n - 2 * Ny + 1, n - Ny + 2, n + Ny + 1};
                    }
                }

                for (int i = 0; i < nn[0].size(); i++)
                {
                    ampo += (0.5 * j2), "S+", nn[0][i], "S-", nn[1][i];
                    ampo += (0.5 * j2), "S-", nn[0][i], "S+", nn[1][i];
                    ampo += (-0.5 * j2), "Nup", nn[0][i], "Ndn", nn[1][i];
                    ampo += (-0.5 * j2), "Ndn", nn[0][i], "Nup", nn[1][i];
                    ampo += -t2, "Cdagup", nn[0][i], "Cup", nn[1][i];
                    ampo += -t2, "Cdagup", nn[1][i], "Cup", nn[0][i];
                    ampo += -t2, "Cdagdn", nn[0][i], "Cdn", nn[1][i];
                    ampo += -t2, "Cdagdn", nn[1][i], "Cdn", nn[0][i];

                    if ((DoMidEnergyPerSite) && ((nn[0][i] - 1) / Ny + 1 > mid_from - 1) && ((nn[0][i] - 1) / Ny + 1 <= mid_to) && ((nn[1][i] - 1) / Ny + 1 > mid_from) && ((nn[1][i] - 1) / Ny + 1 <= mid_to))
                    {
                        ampo_mid += (0.5 * j2), "S+", nn[0][i], "S-", nn[1][i];
                        ampo_mid += (0.5 * j2), "S-", nn[0][i], "S+", nn[1][i];
                        ampo_mid += (-0.5 * j2), "Nup", nn[0][i], "Ndn", nn[1][i];
                        ampo_mid += (-0.5 * j2), "Ndn", nn[0][i], "Nup", nn[1][i];
                        ampo_mid += -t2, "Cdagup", nn[0][i], "Cup", nn[1][i];
                        ampo_mid += -t2, "Cdagup", nn[1][i], "Cup", nn[0][i];
                        ampo_mid += -t2, "Cdagdn", nn[0][i], "Cdn", nn[1][i];
                        ampo_mid += -t2, "Cdagdn", nn[1][i], "Cdn", nn[0][i];
                    }
                }
            }
        }
        if (jk != 0)
        {
            // all clockwise
            //-i/2 S_z S_- S_+
            // i/2 S_z S_+ S_-
            //-i/2 S_+ S_z S_-
            // i/2 S_- S_z S_+
            //-i/2 S_- S_+ S_z
            // i/2 S_+ S_- S_z
            jk = jk / 2.;
            vector<int> bond_120 = {1 - Ny, 1},
                        bond_60 = {1, Ny},
                        bond_120_P = {1 - 2 * Ny, 1 - Ny},
                        bond_60_P = {1 - Ny, Ny};
            for (int n = 1; n <= N; n++)
            {
                int x = (n - 1) / Ny + 1;
                int y = (n - 1) % Ny + 1;
                vector<vector<int>> myvector;
                if (x == 1)
                {
                    if (y == Ny)
                    {
                        if (yperiodic)
                        {
                            myvector.emplace_back(bond_60_P);
                        }
                    }
                    else
                    {
                        myvector.emplace_back(bond_60);
                    }
                }
                else if (x == Nx)
                {
                    if (y == Ny)
                    {
                        if (yperiodic)
                        {
                            myvector.emplace_back(bond_120_P);
                        }
                    }
                    else
                    {
                        myvector.emplace_back(bond_120);
                    }
                }
                else
                {
                    if (y == Ny)
                    {
                        if (yperiodic)
                        {
                            myvector.emplace_back(bond_60_P);
                            myvector.emplace_back(bond_120_P);
                        }
                    }
                    else
                    {
                        myvector.emplace_back(bond_60);
                        myvector.emplace_back(bond_120);
                    }
                }

                for (int i = 0; i < myvector.size(); i++)
                {
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "Sz", n, "S-", n + myvector[i][0], "S+", n + myvector[i][1];
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "Sz", n, "S+", n + myvector[i][0], "S-", n + myvector[i][1];
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "S+", n, "Sz", n + myvector[i][0], "S-", n + myvector[i][1];
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "S-", n, "Sz", n + myvector[i][0], "S+", n + myvector[i][1];
                    ampo += (-(jk)*std::polar(1., M_PI / 2.)), "S-", n, "S+", n + myvector[i][0], "Sz", n + myvector[i][1];
                    ampo += ((jk)*std::polar(1., M_PI / 2.)), "S+", n, "S-", n + myvector[i][0], "Sz", n + myvector[i][1];
                    if ((DoMidEnergyPerSite) && (x > mid_from) && (x <= mid_to))
                    {
                        ampo_mid += (-(jk)*std::polar(1., M_PI / 2.)), "Sz", n, "S-", n + myvector[i][0], "S+", n + myvector[i][1];
                        ampo_mid += ((jk)*std::polar(1., M_PI / 2.)), "Sz", n, "S+", n + myvector[i][0], "S-", n + myvector[i][1];
                        ampo_mid += (-(jk)*std::polar(1., M_PI / 2.)), "S+", n, "Sz", n + myvector[i][0], "S-", n + myvector[i][1];
                        ampo_mid += ((jk)*std::polar(1., M_PI / 2.)), "S-", n, "Sz", n + myvector[i][0], "S+", n + myvector[i][1];
                        ampo_mid += (-(jk)*std::polar(1., M_PI / 2.)), "S-", n, "S+", n + myvector[i][0], "Sz", n + myvector[i][1];
                        ampo_mid += ((jk)*std::polar(1., M_PI / 2.)), "S+", n, "S-", n + myvector[i][0], "Sz", n + myvector[i][1];
                    }
                }
            }
            jk *= 2.;
        }
    }
    else
        Error("Wrong lattice type");
    // edge pining field
    for (int n = 1; n <= Ny; n++)
    {
        int sign = (n % 2 == 1) ? -1 : 1;
        ampo += sign * hz, "Sz", n;
    }
    auto H = toMPO(ampo);
    auto H_mid = toMPO(ampo_mid);
    cout << "MPO Done" << endl;
    // comments on initial MPS
    // string SorC = "";

    cout << "MPS" << endl;
    // Set the initial wavefunction matrix product state
    // to be a Neel state.
    MPS psi(sites);
    auto state = InitState(sites);
    for (int i = 1; i <= N; ++i)
    {
        if (i % 2 == 1)
            state.set(i, "Up");
        else
            state.set(i, "Dn");
    }
    for (int i = 0; i < (doping / 2); ++i)
    {
        state.set(N / 2 - i, "Emp");
        state.set(N / 2 + i + 1, "Emp");
    }

    if (ReadPsiSites)
    {
        readFromFile("psi_tJ_" + lattice_type + "_Nx_" + strNx + "_Ny_" + strNy + "_doping_" + strdoping + "_t1_" + strt1 + "_j2_" + strj2 + "_j2scale_" + strj2_scale + "_jk_" + strjk + CmtIn, psi);
    }
    else
    {
        psi = MPS(state);
        auto psi1 = psi;
        // cout << "isMPScomplexinital = " << isComplex(psi) << endl;
        //
        // complex initial state?
        psi *= 2.0 * Cplx_i;
        psi = sum(psi, psi1, args);
        // creat random initial psi??????
        /*cout << "MPS check 0" << endl;
        psi = applyMPO(H, psi, {"MaxDim=", 50});
        psi.noPrime();
        psi.position(N / 2);
        psi.normalize();
        psi = applyMPO(H, psi, {"MaxDim=", 50});
        psi.noPrime();
        psi.position(N / 2);
        psi.normalize();
        psi = applyMPO(H, psi, {"MaxDim=", 50});
        psi.noPrime();
        psi.position(N / 2);
        psi.normalize();
        for (int j = 1; j <= N; ++j)
        {
            psi.ref(j).randomize({"Complex", true});
        }
        cout << "MPS orthogonal center:" << endl;*/
        psi.position(N / 2);
        cout << "MPS check 2" << endl;
        psi.normalize();
    }
    cout << "initial wave function bond dimension: " << maxLinkDim(psi) << endl;
    cout << "isMPScomplexBeforeDMRG = " << isComplex(psi) << endl;
    println(totalQN(psi));
    cout << "MPS Done" << endl;
    //
    // inner calculates matrix elements of MPO's with respect to MPS's
    // inner(psi,H,psi) = <psi|H|psi>
    //
    printfln("Initial energy = %.8f", innerC(psi, H, psi));

    // total Sz for initial state
    /*auto SZOI = AutoMPO(sites);
    for (int i = 1; i <= N; i++)
    {
        SZOI += "Sz", i;
    }
    auto SZI = toMPO(SZOI);
    cout << "Using AutoMPO initial TotalSz = " << innerC(psi, SZI, psi) << endl;*/
    //
    // Begin the DMRG calculation
    //
    if (DoDMRG)
    {
        cout << "dmrg" << endl;
        auto obs = DMRGObserver(psi, {"EnergyErrgoal", ErrorGoal});
        auto [energy, psi0] = dmrg(H, psi, sweeps, obs, args1);
        psi = psi0;

        // cout << "doing excited states" << endl;
        // auto wfs = std::vector<MPS>(1);
        // wfs.at(0) = psi;
        //
        //  Here the Weight option sets the energy penalty for
        //  psi1 having any overlap with psi
        //
        // auto [en1, psi1] = dmrg(H, wfs, randomMPS(state), sweeps, {"Quiet=", true, "Weight=", 20.0});
        // psi = psi1;
        cout << "dmrg done" << endl;
    }

    // save the ground state to disk
    if (SavePsiSites)
    {
        writeToFile("sites_tJ_" + lattice_type + "_Nx_" + strNx + "_Ny_" + strNy + "_doping_" + strdoping + "_t1_" + strt1 + "_j2_" + strj2 + "_j2scale_" + strj2_scale + "_jk_" + strjk + CmtOut, sites);
        writeToFile("psi_tJ_" + lattice_type + "_Nx_" + strNx + "_Ny_" + strNy + "_doping_" + strdoping + "_t1_" + strt1 + "_j2_" + strj2 + "_j2scale_" + strj2_scale + "_jk_" + strjk + CmtOut, psi);
        cout << "save psi/sites done" << endl;
    }

    if (DoMidDMRGSweeps)
    {
        cout << "Mid dmrg Sweeps from " << mid_from << " to " << mid_to << endl;
        // flip the spin in the center
        psi.position(N / 2);
        auto op_up = sites.op("S+", N / 2);
        psi.set(N / 2, noPrime(psi(N / 2) * op_up));
        auto psi_mid = MPS(((mid_to - mid_from + 1) * Ny));
        auto H_mid1 = MPO(((mid_to - mid_from + 1) * Ny));
        // cout << "length of H_mid1 = " << length(H_mid1) << endl;
        // cout << "length of psi_mid = " << length(psi_mid) << endl;
        //  replace tensor in the middle
        int count4 = 1;
        for (int i = ((mid_from - 1) * Ny + 1); i <= (mid_to * Ny); i++)
        {
            psi_mid.set(count4, psi(i));
            H_mid1.set(count4, H(i));
            count4++;
        }
        // cout << "count4 = " << count4 << endl;
        //  LH and RH
        auto LH = psi(1) * H(1) * dag(prime(psi(1)));
        for (int k = 2; k < ((mid_from - 1) * Ny + 1); ++k)
        {
            LH *= psi(k);
            LH *= H(k);
            LH *= dag(prime(psi(k)));
        }
        // cout << "LH = " << LH << endl;
        auto RH = psi(N) * H(N) * dag(prime(psi(N)));
        for (int k = (N - 1); k > (mid_to * Ny); k--)
        {
            RH *= psi(k);
            RH *= H(k);
            RH *= dag(prime(psi(k)));
        }
        // cout << "RH = " << RH << endl;
        // DMRG sweeps
        auto sweeps_mid = Sweeps({"Nsweep=", 18,
                                  "MaxDim=", 3000,
                                  "MinDim=", 10,
                                  "Noise", 0.,
                                  "Niter", 2,
                                  "Cutoff=", 5E-6});
        auto [energy, psi0] = dmrg(H_mid1, LH, RH, psi_mid, sweeps_mid, {"Quiet", quiet});
        cout << "E1 = " << energy << endl;
        cout << "Mid dmrg Sweeps done" << endl;
    }
    cout << "isMPScomplexAfterDMRG = " << isComplex(psi) << endl;

    // Print the final energy reported by DMRG
    //
    // printfln("\nGround State Energy = %.10f", energy);
    printfln("\nUsing inner = %.10f", innerC(psi, H, psi));
    if (DoMidEnergyPerSite)
    {
        printfln("\nenergy per site on the middle half = %.10f", innerC(psi, H_mid, psi) / ((mid_to - mid_from) * Ny));
    }
    println("\nTotal QN of Ground State = ", totalQN(psi));
    cout << endl;

    int x1, x2, posi;
    double sumCorr, Sab;
    posi = 1;
    psi.position(posi);
    auto C = sites.op("Sz", posi) * psi(posi);
    auto D = sites.op("Sz", posi) * psi(posi);
    auto EE = sites.op("Sz", posi) * psi(posi);
    auto FF = sites.op("Sz", posi) * psi(posi + 1);
    auto GG = sites.op("Sz", posi) * psi(posi);
    auto HH = sites.op("Sz", posi) * psi(posi);
    auto II = sites.op("Sz", posi) * psi(posi);
    auto JJ = sites.op("Sz", posi) * psi(posi);

    if (DoCorrelation)
    {
        // Static correlation function (<S(1,1) * S(L/2,L/2)> - <S(1,1)> * <S(L/2,L/2)>) of electron spin (old)
        cout << "Correlation= <S_i*S_j>" << endl;
        int CorrRange = Nx - 10;  // range
        int CorrStart = Ny * 5 - 2; // first site
        int CorrEnd = Ny * 5;   // last site
        double CorrAve[CorrRange + 1];
        for (int i = 0; i <= (CorrRange); i++)
        {
            CorrAve[i] = 0.0;
        }
        ofstream File;
        File.open("Stat" + TypeofCorrelation + "Corr_tJ_" + lattice_type + "_Nx_" + strNx + "_Ny_" + strNy + "_doping_" + strdoping + "_t1_" + strt1 + "_j2_" + strj2 + "_j2scale_" + strj2_scale + "_jk_" + strjk + CmtOut + ".txt");
        File << "#x ; spin correlation ; hopping correlation (real); hopping correlation (imag); density-density correlation (real)" << endl;
        for (int j = CorrStart; j <= CorrEnd; j++)
        {
            File << "#Correlation i0 at site: " << j << endl;
            for (int i = 1; i <= (CorrRange); i++)
            {
                // correlation function (<S(1,1) * S(L/2,L/2)> - <S(1,1)> * <S(L/2,L/2)>) of electron spin
                x1 = j;
                x2 = x1 + i * Ny;
                Sab = 0.0;
                sumCorr = 0.0;

                // single site magnetization
                /*if ((TypeofCorrelation == "S") || (TypeofCorrelation == "Sz"))
                {
                    posi = x1;
                    psi.position(posi);
                    C = dag(prime(psi(posi), "Site")) * sites.op("Sz", posi) * psi(posi);
                    cout << "<Sz(" << posi << ")> = " << C.cplx().real() << endl;
                    Sab += C.cplx().real();

                    posi = x2;
                    psi.position(posi);
                    C = dag(prime(psi(posi), "Site")) * sites.op("Sz", posi) * psi(posi);
                    //printfln("\nSz = %.10f", C);
                    cout << "<Sz(" << posi << ")> = " << C.cplx().real() << endl;
                    Sab = Sab * C.cplx().real();
                }*/

                // Sz * Sz
                if ((TypeofCorrelation == "S") || (TypeofCorrelation == "Sz"))
                {
                    auto op_ap = sites.op("Nup", x1);
                    auto op_bm = sites.op("Nup", x2);
                    psi.position(x1);
                    auto ir = commonInds(psi(x1), psi(x1 + 1), "Link");
                    C = psi(x1) * op_ap * dag(prime(prime(psi(x1), "Site"), ir));
                    for (int k = x1 + 1; k < x2; ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(x2);
                    C *= op_bm;
                    auto jl = commonInds(psi(x2), psi(x2 - 1), "Link");
                    C *= dag(prime(prime(psi(x2), jl), "Site"));
                    sumCorr += eltC(C).real() / 4.0;

                    op_ap = sites.op("Nup", x1);
                    op_bm = sites.op("Ndn", x2);
                    ir = commonInds(psi(x1), psi(x1 + 1), "Link");
                    C = psi(x1) * op_ap * dag(prime(prime(psi(x1), "Site"), ir));
                    for (int k = x1 + 1; k < x2; ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(x2);
                    C *= op_bm;
                    jl = commonInds(psi(x2), psi(x2 - 1), "Link");
                    C *= dag(prime(prime(psi(x2), jl), "Site"));
                    sumCorr -= eltC(C).real() / 4.0;

                    op_ap = sites.op("Ndn", x1);
                    op_bm = sites.op("Nup", x2);
                    ir = commonInds(psi(x1), psi(x1 + 1), "Link");
                    C = psi(x1) * op_ap * dag(prime(prime(psi(x1), "Site"), ir));
                    for (int k = x1 + 1; k < x2; ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(x2);
                    C *= op_bm;
                    jl = commonInds(psi(x2), psi(x2 - 1), "Link");
                    C *= dag(prime(prime(psi(x2), jl), "Site"));
                    sumCorr -= eltC(C).real() / 4.0;

                    op_ap = sites.op("Ndn", x1);
                    op_bm = sites.op("Ndn", x2);
                    ir = commonInds(psi(x1), psi(x1 + 1), "Link");
                    C = psi(x1) * op_ap * dag(prime(prime(psi(x1), "Site"), ir));
                    for (int k = x1 + 1; k < x2; ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(x2);
                    C *= op_bm;
                    jl = commonInds(psi(x2), psi(x2 - 1), "Link");
                    C *= dag(prime(prime(psi(x2), jl), "Site"));
                    sumCorr += eltC(C).real() / 4.0;
                }

                if ((TypeofCorrelation == "S") || (TypeofCorrelation == "Sxy"))
                {
                    // S+ * S-
                    auto op_ap = sites.op("S+", x1);
                    auto op_bm = sites.op("S-", x2);
                    psi.position(x1);
                    auto ir = commonInds(psi(x1), psi(x1 + 1), "Link");
                    C = psi(x1) * op_ap * dag(prime(prime(psi(x1), "Site"), ir));
                    for (int k = x1 + 1; k < x2; ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(x2);
                    C *= op_bm;
                    auto jl = commonInds(psi(x2), psi(x2 - 1), "Link");
                    C *= dag(prime(prime(psi(x2), jl), "Site"));
                    // cout << "0.5<S+(1) * S-(L/2)> = " << (eltC(C).real() / 2.0) << endl;
                    // cout << "Real and Complex part = " << C.cplx() << " , Cplx should be 0" << endl;
                    //
                    sumCorr += eltC(C).real() / 2.0;

                    // S- * S+
                    auto op_am = sites.op("S-", x1);
                    auto op_bp = sites.op("S+", x2);
                    psi.position(x1);
                    ir = commonInds(psi(x1), psi(x1 + 1), "Link");
                    C = psi(x1) * op_am * dag(prime(prime(psi(x1), "Site"), ir));
                    for (int k = x1 + 1; k < x2; ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(x2);
                    C *= op_bp;
                    jl = commonInds(psi(x2), psi(x2 - 1), "Link");
                    C *= dag(prime(prime(psi(x2), jl), "Site"));
                    // cout << "0.5<S-(1) * S+(L/2)> = " << ( eltC(C).real() / 2.0) << endl;
                    // cout << "Real and Complex part = " << C.cplx() << " , Cplx should be 0" << endl;
                    //
                    sumCorr += eltC(C).real() / 2.0;
                }
                // sumCorr -= Sab;
                cout << "Correlation of [" << TypeofCorrelation << "] (" << x1 << ")(" << x2 << ") = " << sumCorr << endl;
                File << i << " " << sumCorr;
                CorrAve[i] += sumCorr;

                // electron hopping correlation
                //  Cdagup_i * Cup_j = (Adagup_i F_i) F_i+1......F_j-1 (Aup_j)
                psi.position(x1);
                auto ir = commonInds(psi(x1), psi(x1 + 1), "Link");
                C = psi(x1) * sites.op("Adagup*F", x1) * dag(prime(prime(psi(x1), "Site"), ir));
                for (int k = x1 + 1; k < x2; ++k)
                {
                    C *= psi(k);
                    C *= sites.op("F", k);
                    C *= dag(prime(prime(psi(k), "Site"), "Link"));
                }
                C *= psi(x2);
                C *= sites.op("Aup", x2);
                auto jl = commonInds(psi(x2), psi(x2 - 1), "Link");
                C *= dag(prime(prime(psi(x2), jl), "Site"));
                auto sumh = eltC(C);
                // Cdagdn_i * Cdn_j = (Adagdn_i) F_i+1......F_j-1 (F_j Adn_j)
                ir = commonInds(psi(x1), psi(x1 + 1), "Link");
                C = psi(x1) * sites.op("Adagdn", x1) * dag(prime(prime(psi(x1), "Site"), ir));
                for (int k = x1 + 1; k < x2; ++k)
                {
                    C *= psi(k);
                    C *= sites.op("F", k);
                    C *= dag(prime(prime(psi(k), "Site"), "Link"));
                }
                C *= psi(x2);
                C *= sites.op("F*Adn", x2);
                jl = commonInds(psi(x2), psi(x2 - 1), "Link");
                C *= dag(prime(prime(psi(x2), jl), "Site"));
                sumh += eltC(C);
                File << " " << sumh.real() << " " << sumh.imag();

                // density density correlation
                auto op_ap = sites.op("Ntot", x1);
                auto op_bm = sites.op("Ntot", x2);
                psi.position(x1);
                ir = commonInds(psi(x1), psi(x1 + 1), "Link");
                C = psi(x1) * op_ap * dag(prime(prime(psi(x1), "Site"), ir));
                for (int k = x1 + 1; k < x2; ++k)
                {
                    C *= psi(k);
                    C *= dag(prime(psi(k), "Link"));
                }
                C *= psi(x2);
                C *= op_bm;
                jl = commonInds(psi(x2), psi(x2 - 1), "Link");
                C *= dag(prime(prime(psi(x2), jl), "Site"));
                sumCorr = eltC(C).real();
                C = dag(prime(psi(x1), "Site")) * sites.op("Ntot", x1) * psi(x1);
                psi.position(x2);
                C *= dag(prime(psi(x2), "Site")) * sites.op("Ntot", x2) * psi(x2);
                File << " " << sumCorr - eltC(C).real() << endl;
            }
            File << endl;
            File << endl;
            File << endl;
        }
        cout << "Averagr of Correlation of spin: " << endl;
        for (int i = 1; i <= (CorrRange - 1); i++)
        {
            cout << i << " " << CorrAve[i] / (CorrEnd - CorrStart + 1) << endl;
        }
        File.close();
    }

    if (DoSpinStructure)
    {
        ofstream SpinStru;
        ofstream SpinStruDataReal;
        SpinStru.open("SpinStructure_tJ_" + lattice_type + "_Nx_" + strNx + "_Ny_" + strNy + "_doping_" + strdoping + "_t1_" + strt1 + "_j2_" + strj2 + "_j2scale_" + strj2_scale + "_jk_" + strjk + CmtOut + ".txt");
        SpinStruDataReal.open("SpinStruDataRealSpace_tJ_" + lattice_type + "_Nx_" + strNx + "_Ny_" + strNy + "_doping_" + strdoping + "_t1_" + strt1 + "_j2_" + strj2 + "_j2scale_" + strj2_scale + "_jk_" + strjk + CmtOut + ".txt");
        cout << endl;
        int howmuchpoint = Ny * 4;
        // int xx = 84 / Ny; old version
        cout << "New definition: Correlation= <S_i*S_j>" << endl;
        int xx = Nx / 2;
        cout << "range dNx = " << xx << endl;
        int CorrRange = Ny * xx;                                  // range
        int CorrStart = (N - (CorrRange + howmuchpoint)) / 2 + 1; // first site
        // int CorrStart = 1;
        // int CorrEnd = CorrStart + howmuchpoint - 1;			  // last site
        cout << "SpinStructure begin: " << endl;
        cout << "total points = " << howmuchpoint << endl;
        cout << "Correlation Range (dNx * Ny) = " << CorrRange << endl;
        cout << "Correlation Starting point = " << CorrStart << endl;
        cout << endl;
        vector<vector<double>> SpinCorrValue(howmuchpoint);
        for (int i = 0; i < howmuchpoint; i++)
        {
            cout << "Calculating point (" << i << ")" << endl;
            for (int j = 1; j <= CorrRange; j++)
            {

                // correlation function
                x1 = i + CorrStart;
                x2 = j + x1;

                // calculating the spin strucure within a finite size

                // single site magnetization
                /*posi = x1;
                psi.position(posi);
                C = dag(prime(psi(posi), "Site")) * sites.op("Sz", posi) * psi(posi);
                Sab = C.cplx().real();
                posi = x2;
                psi.position(posi);
                C = dag(prime(psi(posi), "Site")) * sites.op("Sz", posi) * psi(posi);
                //printfln("\nSz = %.10f", C);
                Sab = Sab * C.cplx().real();*/
                // Sz * Sz
                auto SzCorr = AutoMPO(sites);
                SzCorr += 1.0, "Nup", x1, "Nup", x2;
                SzCorr += -1.0, "Nup", x1, "Ndn", x2;
                SzCorr += -1.0, "Ndn", x1, "Nup", x2;
                SzCorr += 1.0, "Ndn", x1, "Ndn", x2;
                auto SzCorr1 = toMPO(SzCorr);
                // cout << "<Sz(1) * Sz(L/2)> = " << (C.cplx().real()) << endl;
                // cout << "Real and Complex part = " << C.cplx() << " , Cplx should be 0" << endl;
                //
                sumCorr = innerC(psi, SzCorr1, psi).real();
                // S+ * S-
                auto op_ap = sites.op("S+", x1);
                auto op_bm = sites.op("S-", x2);
                psi.position(x1);
                auto ir = commonInds(psi(x1), psi(x1 + 1), "Link");
                C = psi(x1) * op_ap * dag(prime(prime(psi(x1), "Site"), ir));
                for (int k = x1 + 1; k < x2; ++k)
                {
                    C *= psi(k);
                    C *= dag(prime(psi(k), "Link"));
                }
                C *= psi(x2);
                C *= op_bm;
                auto jl = commonInds(psi(x2), psi(x2 - 1), "Link");
                C *= dag(prime(prime(psi(x2), jl), "Site"));
                // cout << "0.5<S+(1) * S-(L/2)> = " << (C.cplx().real() / 2.0) << endl;
                // cout << "Real and Complex part = " << C.cplx() << " , Cplx should be 0" << endl;
                //
                sumCorr += C.cplx().real() / 2.0;
                // S- * S+
                auto op_am = sites.op("S-", x1);
                auto op_bp = sites.op("S+", x2);
                psi.position(x1);
                ir = commonInds(psi(x1), psi(x1 + 1), "Link");
                C = psi(x1) * op_am * dag(prime(prime(psi(x1), "Site"), ir));
                for (int k = x1 + 1; k < x2; ++k)
                {
                    C *= psi(k);
                    C *= dag(prime(psi(k), "Link"));
                }
                C *= psi(x2);
                C *= op_bp;
                jl = commonInds(psi(x2), psi(x2 - 1), "Link");
                C *= dag(prime(prime(psi(x2), jl), "Site"));
                // cout << "0.5<S-(1) * S+(L/2)> = " << ( C.cplx().real() / 2.0) << endl;
                // cout << "Real and Complex part = " << C.cplx() << " , Cplx should be 0" << endl;
                //
                sumCorr += C.cplx().real() / 2.0;
                // New definition: Correlation= <S_i*S_j>
                // sumCorr -= Sab;
                SpinStruDataReal << j << " " << sumCorr << endl;
                SpinCorrValue[i].emplace_back(sumCorr);
            }
            SpinStruDataReal << endl;
            SpinStruDataReal << endl;
            SpinStruDataReal << endl;
            SpinStruDataReal << endl;
        }
        SpinStruDataReal.close();

        double Sq = 0.0;
        if (lattice_type == "HoneycombXC")
        {
            for (double qx = (-M_PI * 2. / sqrt(3)); qx < (M_PI * 2. / sqrt(3) + 0.05); qx += (4. * M_PI / Ny / sqrt(3)))
            {
                for (double qy = (-4. * M_PI / 3.); qy < (4. * M_PI / 3. + 0.05); qy += (8. * M_PI / Ny / 3.))
                {
                    Sq = 0.0;
                    for (int i = 0; i < howmuchpoint; i++)
                    {
                        for (int j = 0; j < (int)SpinCorrValue[i].size(); j++)
                        {
                            int xx = (j + 1) / Ny;
                            int yy = (j + 1) % Ny;
                            if ((i + CorrStart) % 2 == 0)
                            {
                                if ((i + CorrStart + j + 1) % 2 == 0)
                                {
                                    Sq += 2. * SpinCorrValue[i][j] * cos(qx * (xx * sqrt(3) - (yy * sqrt(3) / 4)) + qy * yy / 2 * 1.5);
                                }
                                else
                                {
                                    Sq += 2. * SpinCorrValue[i][j] * cos(qx * (xx * sqrt(3) - ((yy - 1) * sqrt(3) / 4)) + qy * ((yy - 1) / 2 * 1.5 + 1.));
                                }
                            }
                            else
                            {
                                if ((i + CorrStart + j + 1) % 2 == 0)
                                {
                                    Sq += 2. * SpinCorrValue[i][j] * cos(qx * (xx * sqrt(3) - ((yy + 1) * sqrt(3) / 4)) + qy * ((yy - 1) / 2 * 1.5 + 0.5));
                                }
                                else
                                {
                                    Sq += 2. * SpinCorrValue[i][j] * cos(qx * (xx * sqrt(3) - (yy * sqrt(3) / 4)) + qy * yy / 2 * 1.5);
                                }
                            }
                        }
                    }
                    Sq /= howmuchpoint;
                    SpinStru << qx << " " << qy << " " << Sq << endl;
                }
            }
        }
        else if (lattice_type == "Square")
        {
            for (double qx = (-M_PI * 2.); qx < (M_PI * 2. + 0.02); qx += (2. * M_PI / Ny))
            {
                for (double qy = ((-2.) * M_PI); qy < (2. * M_PI + 0.02); qy += (2. * M_PI / Ny))
                {
                    Sq = 0.0;
                    for (int i = 0; i < howmuchpoint; i++)
                    {
                        for (int j = 0; j < (int)SpinCorrValue[i].size(); j++)
                        {
                            int xx = (j + 1) / Ny;
                            int yy = (j + 1) % Ny;
                            Sq += 2. * SpinCorrValue[i][j] * cos(qx * xx + qy * yy);
                        }
                    }
                    Sq /= howmuchpoint;
                    SpinStru << qx << " " << qy << " " << Sq << endl;
                }
            }
        }
        SpinStru.close();
    }

    //Chirality for Triangle
    if (DoChiralityTriangle)
    {
        double SumValue = 0;
        double SumStagger = 0;
        for (int y = 1; y <= Ny; y++)
        {
            for (int x = (Nx / 4 + 1); x < (Nx / 4 * 3); x++)
            {
                int i = y + Ny * (x - 1);

                if (y < Ny)
                {
                    psi.position(i);
                    //
                    // -i/2 S^z_i S^-_i+1 S^+_i+Ny
                    //
                    auto ir = commonInds(psi(i), psi(i + 1), "Link");
                    C = psi(i) * sites.op("Sz", i) * dag(prime(prime(psi(i), "Site"), ir));
                    C *= psi(i + 1);
                    C *= sites.op("S-", i + 1);
                    C *= dag(prime(psi(i + 1)));
                    for (int k = (i + 2); k < (i + Ny); ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(i + Ny);
                    C *= sites.op("S+", i + Ny);
                    auto jl = commonInds(psi(i + Ny), psi(i + Ny - 1), "Link");
                    C *= dag(prime(prime(psi(i + Ny), jl), "Site"));
                    C = C * (-0.5) * Cplx_i;
                    D = C;
                    //
                    // i/2 S_z S_+ S_-
                    //
                    ir = commonInds(psi(i), psi(i + 1), "Link");
                    C = psi(i) * sites.op("Sz", i) * dag(prime(prime(psi(i), "Site"), ir));
                    C *= psi(i + 1);
                    C *= sites.op("S+", i + 1);
                    C *= dag(prime(psi(i + 1)));
                    for (int k = (i + 2); k < (i + Ny); ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(i + Ny);
                    C *= sites.op("S-", i + Ny);
                    jl = commonInds(psi(i + Ny), psi(i + Ny - 1), "Link");
                    C *= dag(prime(prime(psi(i + Ny), jl), "Site"));
                    C = C * (0.5) * Cplx_i;
                    D += C;
                    //
                    // -i/2 S_+ S_z S_-
                    //
                    ir = commonInds(psi(i), psi(i + 1), "Link");
                    C = psi(i) * sites.op("S+", i) * dag(prime(prime(psi(i), "Site"), ir));
                    C *= psi(i + 1);
                    C *= sites.op("Sz", i + 1);
                    C *= dag(prime(psi(i + 1)));
                    for (int k = (i + 2); k < (i + Ny); ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(i + Ny);
                    C *= sites.op("S-", i + Ny);
                    jl = commonInds(psi(i + Ny), psi(i + Ny - 1), "Link");
                    C *= dag(prime(prime(psi(i + Ny), jl), "Site"));
                    C = C * (-0.5) * Cplx_i;
                    D += C;
                    //
                    // i/2 S_- S_z S_+
                    //
                    ir = commonInds(psi(i), psi(i + 1), "Link");
                    C = psi(i) * sites.op("S-", i) * dag(prime(prime(psi(i), "Site"), ir));
                    C *= psi(i + 1);
                    C *= sites.op("Sz", i + 1);
                    C *= dag(prime(psi(i + 1)));
                    for (int k = (i + 2); k < (i + Ny); ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(i + Ny);
                    C *= sites.op("S+", i + Ny);
                    jl = commonInds(psi(i + Ny), psi(i + Ny - 1), "Link");
                    C *= dag(prime(prime(psi(i + Ny), jl), "Site"));
                    C = C * (0.5) * Cplx_i;
                    D += C;
                    //
                    // -i/2 S_- S_+ S_z
                    //
                    ir = commonInds(psi(i), psi(i + 1), "Link");
                    C = psi(i) * sites.op("S-", i) * dag(prime(prime(psi(i), "Site"), ir));
                    C *= psi(i + 1);
                    C *= sites.op("S+", i + 1);
                    C *= dag(prime(psi(i + 1)));
                    for (int k = (i + 2); k < (i + Ny); ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(i + Ny);
                    C *= sites.op("Sz", i + Ny);
                    jl = commonInds(psi(i + Ny), psi(i + Ny - 1), "Link");
                    C *= dag(prime(prime(psi(i + Ny), jl), "Site"));
                    C = C * (-0.5) * Cplx_i;
                    D += C;
                    //
                    // i/2 S_+ S_- S_z
                    //
                    ir = commonInds(psi(i), psi(i + 1), "Link");
                    C = psi(i) * sites.op("S+", i) * dag(prime(prime(psi(i), "Site"), ir));
                    C *= psi(i + 1);
                    C *= sites.op("S-", i + 1);
                    C *= dag(prime(psi(i + 1)));
                    for (int k = (i + 2); k < (i + Ny); ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(i + Ny);
                    C *= sites.op("Sz", i + Ny);
                    jl = commonInds(psi(i + Ny), psi(i + Ny - 1), "Link");
                    C *= dag(prime(prime(psi(i + Ny), jl), "Site"));
                    C = C * (0.5) * Cplx_i;
                    D += C;
                    cout << D.cplx().real() << " ";
                }
                else
                {
                    psi.position(i - Ny + 1);
                    //
                    // -i/2 S^z_i S^-_i-Ny+1 S^+_i+Ny
                    //
                    auto ir = commonInds(psi(i - Ny + 1), psi(i - Ny + 2), "Link");
                    C = psi(i - Ny + 1) * sites.op("S-", i - Ny + 1) * dag(prime(prime(psi(i - Ny + 1), "Site"), ir));
                    for (int k = (i - Ny + 2); k < i; ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(i);
                    C *= sites.op("Sz", i);
                    C *= dag(prime(psi(i)));
                    for (int k = (i + 1); k < (i + Ny); ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(i + Ny);
                    C *= sites.op("S+", i + Ny);
                    auto jl = commonInds(psi(i + Ny), psi(i + Ny - 1), "Link");
                    C *= dag(prime(prime(psi(i + Ny), jl), "Site"));
                    C = C * (-0.5) * Cplx_i;
                    D = C;
                    //
                    // i/2 S_z S_+ S_-
                    //
                    ir = commonInds(psi(i - Ny + 1), psi(i - Ny + 2), "Link");
                    C = psi(i - Ny + 1) * sites.op("S+", i - Ny + 1) * dag(prime(prime(psi(i - Ny + 1), "Site"), ir));
                    for (int k = (i - Ny + 2); k < i; ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(i);
                    C *= sites.op("Sz", i);
                    C *= dag(prime(psi(i)));
                    for (int k = (i + 1); k < (i + Ny); ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(i + Ny);
                    C *= sites.op("S-", i + Ny);
                    jl = commonInds(psi(i + Ny), psi(i + Ny - 1), "Link");
                    C *= dag(prime(prime(psi(i + Ny), jl), "Site"));
                    C = C * (0.5) * Cplx_i;
                    D += C;
                    //
                    // -i/2 S_+ S_z S_-
                    //
                    ir = commonInds(psi(i - Ny + 1), psi(i - Ny + 2), "Link");
                    C = psi(i - Ny + 1) * sites.op("Sz", i - Ny + 1) * dag(prime(prime(psi(i - Ny + 1), "Site"), ir));
                    for (int k = (i - Ny + 2); k < i; ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(i);
                    C *= sites.op("S+", i);
                    C *= dag(prime(psi(i)));
                    for (int k = (i + 1); k < (i + Ny); ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(i + Ny);
                    C *= sites.op("S-", i + Ny);
                    jl = commonInds(psi(i + Ny), psi(i + Ny - 1), "Link");
                    C *= dag(prime(prime(psi(i + Ny), jl), "Site"));
                    C = C * (-0.5) * Cplx_i;
                    D += C;
                    //
                    // i/2 S_- S_z S_+
                    //
                    ir = commonInds(psi(i - Ny + 1), psi(i - Ny + 2), "Link");
                    C = psi(i - Ny + 1) * sites.op("Sz", i - Ny + 1) * dag(prime(prime(psi(i - Ny + 1), "Site"), ir));
                    for (int k = (i - Ny + 2); k < i; ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(i);
                    C *= sites.op("S-", i);
                    C *= dag(prime(psi(i)));
                    for (int k = (i + 1); k < (i + Ny); ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(i + Ny);
                    C *= sites.op("S+", i + Ny);
                    jl = commonInds(psi(i + Ny), psi(i + Ny - 1), "Link");
                    C *= dag(prime(prime(psi(i + Ny), jl), "Site"));
                    C = C * (0.5) * Cplx_i;
                    D += C;
                    //
                    // -i/2 S_- S_+ S_z
                    //
                    ir = commonInds(psi(i - Ny + 1), psi(i - Ny + 2), "Link");
                    C = psi(i - Ny + 1) * sites.op("S+", i - Ny + 1) * dag(prime(prime(psi(i - Ny + 1), "Site"), ir));
                    for (int k = (i - Ny + 2); k < i; ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(i);
                    C *= sites.op("S-", i);
                    C *= dag(prime(psi(i)));
                    for (int k = (i + 1); k < (i + Ny); ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(i + Ny);
                    C *= sites.op("Sz", i + Ny);
                    jl = commonInds(psi(i + Ny), psi(i + Ny - 1), "Link");
                    C *= dag(prime(prime(psi(i + Ny), jl), "Site"));
                    C = C * (-0.5) * Cplx_i;
                    D += C;
                    //
                    // i/2 S_+ S_- S_z
                    //
                    ir = commonInds(psi(i - Ny + 1), psi(i - Ny + 2), "Link");
                    C = psi(i - Ny + 1) * sites.op("S-", i - Ny + 1) * dag(prime(prime(psi(i - Ny + 1), "Site"), ir));
                    for (int k = (i - Ny + 2); k < i; ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(i);
                    C *= sites.op("S+", i);
                    C *= dag(prime(psi(i)));
                    for (int k = (i + 1); k < (i + Ny); ++k)
                    {
                        C *= psi(k);
                        C *= dag(prime(psi(k), "Link"));
                    }
                    C *= psi(i + Ny);
                    C *= sites.op("Sz", i + Ny);
                    jl = commonInds(psi(i + Ny), psi(i + Ny - 1), "Link");
                    C *= dag(prime(prime(psi(i + Ny), jl), "Site"));
                    C = C * (0.5) * Cplx_i;
                    D += C;
                    cout << D.cplx().real() << " ";
                }
                SumValue += D.cplx().real();
                SumStagger += (y % 2 == 1) ? -(D.cplx().real()) : D.cplx().real();
            }
            cout << endl;
        }
        cout << endl;
        cout << "Sum of middle half Chiral order (+45) = " << (SumValue / (Nx / 2. * Ny)) << endl;
        cout << "Stagger Sum of middle half Chiral order (+45) = " << (SumStagger / (Nx / 2. * Ny)) << endl;
    }

    if (SCcorr)
    {
        ofstream SCf;
        cout << "SC correlation begin" << endl;
        SCf.open("SCCorrelation_tJ_" + lattice_type + "_Nx_" + strNx + "_Ny_" + strNy + "_doping_" + strdoping + "_t1_" + strt1 + "_j2_" + strj2 + "_j2scale_" + strj2_scale + "_jk_" + strjk + CmtOut + ".txt");
        SCf << "#x, bb, ba, bc, aa, ac, cc (singlet pairing)" << endl;
        x1 = Ny * 5 - 2;
        SCf << "#Correlation i0 at site: " << x1 << endl;
        for (int j = 2; j <= (Nx - 10); j++)
        {
            // SC correlation function
            x2 = x1 + j * Ny;
            ///////////////////////////////////////////////////////////////
            // x, a-bond value, a-bond phase, b-bond value, b-bond phase, c-bond value, c-bond phase (singlet pairing)
            vector<int> types = {1, Ny, 1 - Ny};
            SCf << j;
            for (int i = 0; i < types.size(); i++)
            {
                for (int b = i; b < types.size(); b++)
                {
                    // Cdagup_x1 Cdagdn_x1 + Cdagdn_x1 Cdagup_x1
                    auto SCterm = AutoMPO(sites);
                    if (types[i] == 0)
                    {
                        SCterm += 1.0, "Cdagup*Cdagdn", x1, "Cdn*Cup", x2;
                    }
                    else
                    {
                        SCterm += 0.5, "Cdagup", x1, "Cdagdn", x1 + types[i], "Cdn", x2 + types[b], "Cup", x2;
                        SCterm += -0.5, "Cdagup", x1, "Cdagdn", x1 + types[i], "Cup", x2 + types[b], "Cdn", x2;
                        SCterm += -0.5, "Cdagdn", x1, "Cdagup", x1 + types[i], "Cdn", x2 + types[b], "Cup", x2;
                        SCterm += 0.5, "Cdagdn", x1, "Cdagup", x1 + types[i], "Cup", x2 + types[b], "Cdn", x2;
                        // SCterm += 1.0, "Cdagup", x2, "Cdagdn", x2 + types[b], "Cdn", x1 + types[i], "Cup", x1;
                    }
                    auto SCterm1 = toMPO(SCterm);
                    // Cdagup_i * Cdn_j+1 = (Adagup_i F_i) F_i+1......F_j (F_j+1 Adn_j+1)
                    /*psi.position(x1);
                auto ir = commonInds(psi(x1), psi(x1 + 1), "Link");
                C = psi(x1) * sites.op("Adagup*F", x1) * dag(prime(prime(psi(x1), "Site"), ir));
                for (int k = x1 + 1; k < (x2 + types[i]); ++k)
                {
                    C *= psi(k);
                    C *= sites.op("F", k);
                    C *= dag(prime(prime(psi(k), "Link"), "Site"));
                }
                C *= psi(x2 + types[i]);
                C *= sites.op("F*Adn", x2 + types[i]);
                auto jl = commonInds(psi(x2 + types[i]), psi(x2 + types[i] - 1), "Link");
                C *= dag(prime(prime(psi(x2 + types[i]), jl), "Site"));
                auto Sab = C.cplx();
                if (i == 2)
                {
                    cout << "Cdagup Cdn = " << Sab << endl;
                }
                // Cdagdn_i+1 * Cup_j = (Adagdn_i+1) F_i+2......F_j-1 (Aup_j) Adagdn Aup
                psi.position(x1 + types[i]);
                ir = commonInds(psi(x1 + types[i]), psi(x1 + types[i] + 1), "Link");
                C = psi(x1 + types[i]) * sites.op("Adagdn", x1 + types[i]) * dag(prime(prime(psi(x1 + types[i]), "Site"), ir));
                for (int k = x1 + types[i] + 1; k < (x2); ++k)
                {
                    C *= psi(k);
                    C *= sites.op("F", k);
                    C *= dag(prime(prime(psi(k), "Link"), "Site"));
                }
                C *= psi(x2);
                C *= sites.op("Aup", x2);
                jl = commonInds(psi(x2), psi(x2 - 1), "Link");
                C *= dag(prime(prime(psi(x2), jl), "Site"));
                Sab *= C.cplx();*/
                    // Cdagup_i * Cup_j = (Adagup_i F_i) F_i+1......F_j-1 (Aup_j)
                    /*psi.position(x1);
                    auto ir = commonInds(psi(x1), psi(x1 + 1), "Link");
                    C = psi(x1) * sites.op("Adagup*F", x1) * dag(prime(prime(psi(x1), "Site"), ir));
                    for (int k = x1 + 1; k < (x2); ++k)
                    {
                        C *= psi(k);
                        C *= sites.op("F", k);
                        C *= dag(prime(prime(psi(k), "Link"), "Site"));
                    }
                    C *= psi(x2);
                    C *= sites.op("Aup", x2);
                    auto jl = commonInds(psi(x2), psi(x2 - 1), "Link");
                    C *= dag(prime(prime(psi(x2), jl), "Site"));
                    auto Scd = eltC(C);
                    // Cdagdn_i+type[i] * Cdn_j+type[b] = (Adagdn_i+1) F_i+2......F_j (F_i+1 Adn_j+1)
                    psi.position(x1 + types[i]);
                    ir = commonInds(psi(x1 + types[i]), psi(x1 + types[i] + 1), "Link");
                    C = psi(x1 + types[i]) * sites.op("Adagdn", x1 + types[i]) * dag(prime(prime(psi(x1 + types[i]), "Site"), ir));
                    for (int k = x1 + types[i] + 1; k < (x2 + types[b]); ++k)
                    {
                        C *= psi(k);
                        C *= sites.op("F", k);
                        C *= dag(prime(prime(psi(k), "Link"), "Site"));
                    }
                    C *= psi(x2 + types[b]);
                    C *= sites.op("F*Adn", x2 + types[b]);
                    jl = commonInds(psi(x2 + types[b]), psi(x2 + types[b] - 1), "Link");
                    C *= dag(prime(prime(psi(x2 + types[b]), jl), "Site"));
                    Scd *= eltC(C);*/
                    // Cdn_i+1 * Cdagup_j = (Adn_i+1) F_i+2......F_j-1 (Adagup_j)
                    /*psi.position(x1 + types[i]);
                ir = commonInds(psi(x1 + types[i]), psi(x1 + types[i] + 1), "Link");
                C = psi(x1 + types[i]) * sites.op("Adn", x1 + types[i]) * dag(prime(prime(psi(x1 + types[i]), "Site"), ir));
                for (int k = x1 + types[i] + 1; k < (x2); ++k)
                {
                    C *= psi(k);
                    C *= sites.op("F", k);
                    C *= dag(prime(prime(psi(k), "Link"), "Site"));
                }
                C *= psi(x2);
                C *= sites.op("Adagup", x2);
                jl = commonInds(psi(x2), psi(x2 - 1), "Link");
                C *= dag(prime(prime(psi(x2), jl), "Site"));
                auto Sef = C.cplx();
                // Cup_i * Cdagdn_j+1 = (Aup_i F_i) F_i+1......F_j (F_i+1 Adagdn_j+1)
                psi.position(x1);
                ir = commonInds(psi(x1), psi(x1 + 1), "Link");
                C = psi(x1) * sites.op("Aup*F", x1) * dag(prime(prime(psi(x1), "Site"), ir));
                for (int k = x1 + 1; k < (x2 + types[i]); ++k)
                {
                    C *= psi(k);
                    C *= sites.op("F", k);
                    C *= dag(prime(prime(psi(k), "Link"), "Site"));
                }
                C *= psi(x2 + types[i]);
                C *= sites.op("F*Adagdn", x2 + types[i]);
                jl = commonInds(psi(x2 + types[i]), psi(x2 + types[i] - 1), "Link");
                C *= dag(prime(prime(psi(x2 + types[i]), jl), "Site"));
                Sef *= C.cplx();*/

                    // Cup_i * Cdagup_j = (Aup_i F_i) F_i+1......F_j-1 (Adagup_j)
                    /*psi.position(x1);
                    ir = commonInds(psi(x1), psi(x1 + 1), "Link");
                    C = psi(x1) * sites.op("Aup*F", x1) * dag(prime(prime(psi(x1), "Site"), ir));
                    for (int k = x1 + 1; k < (x2); ++k)
                    {
                        C *= psi(k);
                        C *= sites.op("F", k);
                        C *= dag(prime(prime(psi(k), "Link"), "Site"));
                    }
                    C *= psi(x2);
                    C *= sites.op("Adagup", x2);
                    jl = commonInds(psi(x2), psi(x2 - 1), "Link");
                    C *= dag(prime(prime(psi(x2), jl), "Site"));
                    auto Sgh = eltC(C);
                    // Cdn_i+1 * Cdagdn_j+1 = (Adn_i+1) F_i+2......F_j (F_i+1 Adagdn_j+1)
                    psi.position(x1 + types[i]);
                    ir = commonInds(psi(x1 + types[i]), psi(x1 + types[i] + 1), "Link");
                    C = psi(x1 + types[i]) * sites.op("Adn", x1 + types[i]) * dag(prime(prime(psi(x1 + types[i]), "Site"), ir));
                    for (int k = x1 + types[i] + 1; k < (x2 + types[b]); ++k)
                    {
                        C *= psi(k);
                        C *= sites.op("F", k);
                        C *= dag(prime(prime(psi(k), "Link"), "Site"));
                    }
                    C *= psi(x2 + types[b]);
                    C *= sites.op("F*Adagdn", x2 + types[b]);
                    jl = commonInds(psi(x2 + types[b]), psi(x2 + types[b] - 1), "Link");
                    C *= dag(prime(prime(psi(x2 + types[b]), jl), "Site"));
                    Sgh *= eltC(C);*/
                    SCf << " " << abs(innerC(psi, SCterm1, psi)) << " " << arg(innerC(psi, SCterm1, psi));
                }
            }
            SCf << endl;
        }
        SCf << endl;
        SCf.close();
        cout << "SC correlation done" << endl;
    }

    if (DoLocalSpinValue)
    {
        // sttagered magnetization
        // total magnetization
        //
        double sumSz = 0.0;
        double SttaMag = 0.0;
        double ColliMag = 0.0;
        cout << endl;
        cout << "Site position; charge" << endl;
        for (int i = 1; i <= N; i++)
        {
            psi.position(i);
            C = dag(prime(psi(i), "Site")) * sites.op("Ntot", i) * psi(i);
            cout << i << " " << eltC(C).real() << endl;
        }
        cout << endl;
        cout << "Site position; spin" << endl;
        for (int i = 1; i <= N; i++)
        {
            psi.position(i);
            C = dag(prime(psi(i), "Site")) * sites.op("Nup", i) * psi(i);
            auto C1 = dag(prime(psi(i), "Site")) * sites.op("Ndn", i) * psi(i);
            cout << i << " " << eltC(C - C1).real() / 2.0 << endl;
            double sign = ((((i - 1) / Ny) % 2) == 1) ? 1 : -1;
            double Mid = ((((i - 1) / Ny + 1) > (Nx / 4)) && (((i - 1) / Ny + 1) <= (Nx - (Nx / 4)))) ? 1 : 0;
            SttaMag += ((i % 2) == 1) ? (sign * Mid * eltC(C).real()) : ((-sign) * Mid * eltC(C).real());
            ColliMag += ((i % 2) == 1) ? (Mid * eltC(C).real()) : ((-Mid) * eltC(C).real());
            sumSz += eltC(C).real();
        }
        cout << endl;
        cout << "total magnetization = " << sumSz << endl;
        cout << "total magnetization per site = " << (sumSz / N) << endl;
        cout << "Mid half " << (Nx - (Nx / 4) * 2) << " sites in X direction" << endl;
        cout << "staggered magnetization per site (Mid half) = " << (SttaMag / (Ny * (Nx - (Nx / 4) * 2))) << endl;
        cout << "Collinear magnetization per site (Mid half) = " << (ColliMag / (Ny * (Nx - (Nx / 4) * 2))) << endl;
        // spin value
        cout << endl;
        cout << "spin value in x direction (y=Ny/2): " << endl;
        int count = 1;
        for (int i = (Ny / 2); i <= N; i = i + Ny)
        {
            psi.position(i);
            C = dag(prime(psi(i), "Site")) * sites.op("Nup", i) * psi(i);
            auto C1 = dag(prime(psi(i), "Site")) * sites.op("Ndn", i) * psi(i);
            cout << count << " " << eltC(C - C1).real() / 2.0 << endl;
            count++;
        }
    }
    // VEE calculation
    if (DoVEE)
    {
        cout << endl;
        /*cout << " site ; VEE (revised)" << endl;
        //int b = N / 2;
        for (int b = 1; b < N; b++)
        {
            psi.position(b);
            //SVD this wavefunction to get the spectrum
            //of density-matrix eigenvalues
            ITensor wf = psi(b) * psi(b + 1);
            auto U = psi(b);
            ITensor S, V;
            auto spectrum = svd(wf, U, S, V, {"ComputeQNs", true});

            //Apply von Neumann formula
            //to the squares of the singular values
            //Real SvN = 0.;
            //for (auto eig : spectrum.eigsKept())
            //{
            //    if (eig > 1E-12)
            //        SvN += -eig * log(eig);
            //}

            auto l = leftLinkIndex(psi, b);
            auto s = siteIndex(psi, b);
            auto [U, S, V] = svd(psi(b), {l, s});
            auto u = commonIndex(U, S);

            //Apply von Neumann formula
            //to the squares of the singular values
            Real SvN = 0.;
            for (auto n : range1(dim(u)))
            {
                auto Sn = elt(S, n, n);
                auto p = sqr(Sn);
                if (p > 1E-12)
                    SvN += -p * log(p);
            }

            cout << b << " " << SvN << endl;
        }*/
        cout << endl;
        cout << " spectrum in the middle:" << endl;
        int x = N / 2;
        psi.position(x);
        ITensor wf = psi(x) * psi(x + 1);
        auto U = psi(x);
        ITensor S, V;
        auto spectrum = svd(wf, U, S, V, {"ComputeQNs", true});

        // Apply von Neumann formula
        // to the squares of the singular values
        Real SvN = 0.;
        for (int i = 0; i < spectrum.eigsKept().size(); i++)
        {
            if (spectrum.eigsKept()[i] > 1E-6)
            {
                cout << spectrum.qns()[i].val("Nf") << " " << spectrum.qns()[i].val("Sz") << " " << spectrum.eigsKept()[i] << endl;
                SvN += -spectrum.eigsKept()[i] * log(spectrum.eigsKept()[i]);
            }
        }
        cout << endl;
        cout << endl;
        cout << " site ; VEE (revised)" << endl;
        cout << x << " " << SvN << endl;

        // finding the largest eigenvalue with another way for checking
        /*double first = 0, second = 0, third = 0, fourth = 0;
        for (auto p : spectrum.eigs())
        {
            if (p > first)
                first = p;
        }
        //finding the second largest eigenvalue
        for (auto p : spectrum.eigs())
        {
            if (p < first)
            {
                if (p > second)
                    second = p;
            }
        }
        //finding the third largest eigenvalue
        for (auto p : spectrum.eigs())
        {
            if ((p < first) && (p < second))
            {
                if (p > third)
                    third = p;
            }
        }
        //finding the fourth largest eigenvalue
        for (auto p : spectrum.eigs())
        {
            if ((p < first) && (p < second) && (p < third))
            {
                if (p > fourth)
                    fourth = p;
            }
        }
        cout << "traditional way:" << endl;
        cout << first << " " << second << " " << third << " " << fourth << endl;*/
    }

    // Entanglement Spectrum calculation

    // print parameters
    //
    cout << " j1 = " << j1 << endl;
    cout << " j2 = " << j2 << endl;
    cout << " j2_scale = " << j2_scale << endl;
    cout << " jk = " << jk << endl;
    cout << " t1 = " << t1 << endl;
    cout << " doping number = " << doping << endl;
    cout << " Nx = " << Nx << endl;
    cout << " Ny = " << Ny << endl;
    cout << " hz = " << hz << endl;
    cout << " jkedge = " << jkedge << endl;
    cout << " lattice type = " << lattice_type << endl;

    /*psi.position(N / 2);
    for (int i = 1; i < 4; i++)
    {
        cout << i << "th Tensor " << endl;
        print(psi(i));
        cout << N - i + 1 << "th Tensor " << endl;
        print(psi(N - i + 1));
    }*/

    // double E0 = inner(psi, H, psi);
    // double E02 = inner(psi, H, H, psi);
    // cout << " Ground State variance = " << E02 - E0 * E0 << endl;
    cout << endl;
    return 0;
}
