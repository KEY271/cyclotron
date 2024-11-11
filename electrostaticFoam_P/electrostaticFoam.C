#include "fvCFD.H"
#include "interpolationCellPoint.H"

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for electrostatics."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "\nInitializing\n" << endl;

    const volVectorField E1 = (-fvc::grad(phi1) * V).ref();
    const volVectorField E2 = (-fvc::grad(phi2) * Vf).ref();
    const volVectorField E3 = (-fvc::grad(phi3) * Vi / 2.0).ref();
    const volVectorField E23 = (E2 + E3).ref();
    const interpolationCellPoint<vector> E1Interp(E1);
    const interpolationCellPoint<vector> E23Interp(E23);

    const double timestep = 0.001e-6;
    double time = 0;

    const fileName dir = "results";
    if (!isDir(dir)) {
        mkDir(dir);
    }

    Random rng;

    Info << "\nCalculating\n" << endl;

    const fileName file = dir / "result.csv";
    OFstream os(file);

    for (double b = 0; b < 0.3; b += 0.005) {
        int count = 0;
        for (int n = 1; n <= 1000; ++n) {
            // const fileName file = dir / "trajectory" + std::to_string(n);
            // OFstream os(file);
            scalar angle = 2.0 * M_PI * rng.sample01<scalar>();
            scalar r = 0.0005 + 0.003 * rng.sample01<scalar>();
            scalar theta = rng.sample01<scalar>();
            point p(0.005 + r * std::cos(angle), -0.005 + r * std::sin(angle), 0.01 + 0.004 * 2.0 * (rng.sample01<scalar>() - 0.5));
            vector v(0, 0, 0);
            for (int i = 0; i < 100000; ++i) {
                if ((p.x() - 0.005) * (p.x() - 0.005) + (p.y() + 0.005) * (p.y() + 0.005) < 0.00025 * 0.00025 &&
                    p.z() > 0.006 && p.z() < 0.014) {
                    break;
                }
                if (p.x() * p.x() + p.y() * p.y() > 0.025 * 0.025) {
                    break;
                }
                if (p.z() < 0.004 || p.z() > 0.014) {
                    break;
                }
                if (p.z() > 0.006 && p.z() < 0.014 && p.x() > -0.002 && p.x() < 0.002 && p.y() < -0.015 && p.y() > -0.025) {
                    count++;
                    break;
                }
                label cell = mesh.findCell(p);
                if (cell == -1) {
                    cell = mesh.findNearestCell(p);
                }
                const vector dv1 = (
                    E1Interp.interpolate(p, cell) * std::sin(2.0 * M_PI * (f.value() * time + theta))
                    + E23Interp.interpolate(p, cell)
                    + vector(v.y() * b, -v.x() * b, 0.0)) * (e / m * timestep).value();

                const point p2 = p + 0.5 * (v + dv1) * timestep;
                const vector dv2 = (
                    E1Interp.interpolate(p2, cell) * std::sin(2.0 * M_PI * (f.value() * (time + timestep / 2.0) + theta))
                    + E23Interp.interpolate(p2, cell)
                    + vector(v.y() * b, -v.x() * b, 0.0)) * (e / m * timestep).value();

                const point p3 = p2 + 0.5 * (v + dv2) * timestep;
                const vector dv3 = (
                    E1Interp.interpolate(p3, cell) * std::sin(2.0 * M_PI * (f.value() * (time + timestep / 2.0) + theta))
                    + E23Interp.interpolate(p3, cell)
                    + vector(v.y() * b, -v.x() * b, 0.0)) * (e / m * timestep).value();

                const point p4 = p3 + (v + dv3) * timestep;
                const vector dv4 = (
                    E1Interp.interpolate(p4, cell) * std::sin(2.0 * M_PI * (f.value() * (time + timestep) + theta))
                    + E23Interp.interpolate(p4, cell)
                    + vector(v.y() * b, -v.x() * b, 0.0)) * (e / m * timestep).value();

                p += (v + (dv1 * 2.0 + dv2 * 2.0 + dv3) / 6.0) * timestep;
                v += (dv1 + dv2 * 2.0 + dv3 * 2.0 + dv4) / 6.0;

                // os << time << " " << p.x() << " " << p.y() << " " << p.z() << " " << (v & v * m / 2 / e).value() << endl;
                time += timestep;
            }
            time = 0;
        }
        os << b << " " << count << endl;
        Info << "B " << b << " T done; count = " << count << "\n" << endl;
    }

    Info << "End\n" << endl;

    return 0;
}
