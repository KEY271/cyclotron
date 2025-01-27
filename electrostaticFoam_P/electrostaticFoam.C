#include "fvCFD.H"
#include "interpolationCellPoint.H"

int main(int argc, char *argv[]) {
    argList::addNote(
        "Solver for electrostatics."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "Initializing" << endl;

    const volVectorField E1 = (-fvc::grad(phi1) * V).ref();
    const volVectorField E2 = (-fvc::grad(phi2) * Vf).ref();
    const volVectorField E3 = (-fvc::grad(phi3) * Vi / 2.0).ref();
    const volVectorField E4 = (-fvc::grad(phi4) * Vc).ref();
    const interpolationCellPoint<vector> E1Interp(E1);
    const interpolationCellPoint<vector> E2Interp(E2);
    const interpolationCellPoint<vector> E3Interp(E3);
    const interpolationCellPoint<vector> E4Interp(E4);

    const double timestep = 0.001e-6;
    double time = 0;

    const fileName dir = "results";
    if (!isDir(dir)) {
        mkDir(dir);
    }
    const fileName file = dir / "result";
    OFstream os(file);
    os << "b count fallen r_hit z_hit" << endl;

    Random rng;

    Info << "Starting" << endl;

    for (double b = 0; b < 0.301; b += 0.01) {
        int count = 0;
        int fallen = 0;
        int r_hit = 0;
        int z_hit = 0;
        for (int n = 1; n <= 1000; ++n) {
            scalar angle = 2.0 * M_PI * rng.sample01<scalar>();
            scalar r = 0.0005 + 0.003 * rng.sample01<scalar>();
            scalar theta = rng.sample01<scalar>();
            point p(r * std::cos(angle), -0.005 + r * std::sin(angle), 0.01 + 0.004 * 2.0 * (rng.sample01<scalar>() - 0.5));
            vector v(0, 0, 0);
            for (int i = 0; i < 100000; ++i) {
                if (p.x() * p.x() + (p.y() + 0.005) * (p.y() + 0.005) < 0.00025 * 0.00025 &&
                    p.z() > 0.006 && p.z() < 0.014) {
                    fallen++;
                    break;
                }
                if (p.x() * p.x() + (p.z() - 0.006) * (p.z() - 0.006) < 0.00025 * 0.00025 &&
                    p.y() > -0.004) {
                    fallen++;
                    break;
                }
                if (p.x() * p.x() + (p.z() - 0.014) * (p.z() - 0.014) < 0.00025 * 0.00025 &&
                    p.y() > -0.004) {
                    fallen++;
                    break;
                }
                if ((p.x() + 0.005) * (p.x() + 0.005) + p.y() * p.y() > 0.025 * 0.025 && p.x() < -0.005) {
                    r_hit++;
                    break;
                }
                if (p.x() * p.x() + p.y() * p.y() > 0.037 * 0.037) {
                    r_hit++;
                    break;
                }
                if (p.z() < 0.004 || p.z() > 0.014) {
                    z_hit++;
                    break;
                }
                if (p.z() > 0.006 && p.z() < 0.014 && p.y() > -0.002 && p.y() < 0.002 && p.x() > 0.018 && p.x() < 0.032) {
                    count++;
                    break;
                }
                label cell = mesh.findCell(p);
                if (cell == -1) {
                    cell = mesh.findNearestCell(p);
                }
                const vector dv1 = (
                    E1Interp.interpolate(p, cell) * std::sin(2.0 * M_PI * (f.value() * time + theta))
                    + E2Interp.interpolate(p, cell)
                    + E3Interp.interpolate(p, cell)
                    + E4Interp.interpolate(p, cell)
                    + vector(v.y() * b, -v.x() * b, 0.0)) * (e / m * timestep).value();

                const point p2 = p + 0.5 * (v + dv1) * timestep;
                const vector dv2 = (
                    E1Interp.interpolate(p2, cell) * std::sin(2.0 * M_PI * (f.value() * (time + timestep / 2.0) + theta))
                    + E2Interp.interpolate(p2, cell)
                    + E3Interp.interpolate(p2, cell)
                    + E4Interp.interpolate(p2, cell)
                    + vector(v.y() * b, -v.x() * b, 0.0)) * (e / m * timestep).value();

                const point p3 = p2 + 0.5 * (v + dv2) * timestep;
                const vector dv3 = (
                    E1Interp.interpolate(p3, cell) * std::sin(2.0 * M_PI * (f.value() * (time + timestep / 2.0) + theta))
                    + E2Interp.interpolate(p3, cell)
                    + E3Interp.interpolate(p3, cell)
                    + E4Interp.interpolate(p3, cell)
                    + vector(v.y() * b, -v.x() * b, 0.0)) * (e / m * timestep).value();

                const point p4 = p3 + (v + dv3) * timestep;
                const vector dv4 = (
                    E1Interp.interpolate(p4, cell) * std::sin(2.0 * M_PI * (f.value() * (time + timestep) + theta))
                    + E2Interp.interpolate(p4, cell)
                    + E3Interp.interpolate(p4, cell)
                    + E4Interp.interpolate(p4, cell)
                    + vector(v.y() * b, -v.x() * b, 0.0)) * (e / m * timestep).value();

                p += (v + (dv1 * 2.0 + dv2 * 2.0 + dv3) / 6.0) * timestep;
                v += (dv1 + dv2 * 2.0 + dv3 * 2.0 + dv4) / 6.0;
                time += timestep;
            }
            time = 0;
            Info << n << "\r" << flush;
        }
        Info << endl;

        os << b << " " << count << " " << fallen << " " << r_hit << " " << z_hit << endl;
        os.flush();
        Info << "b = " << b << ", count = " << count << ", fallen = " << fallen << ", r_hit = " << r_hit << ", z_hit = " << z_hit << endl;
    }

    Info << "Finished" << endl;

    return 0;
}
