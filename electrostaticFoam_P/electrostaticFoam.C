#include "fvCFD.H"
#include "interpolationCellPoint.H"

int main(int argc, char *argv[]) {
    argList::addNote(
        "Solver for charged particles."
    );

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info << "Initializing" << endl;

    // 電場
    const volVectorField E1 = (-fvc::grad(phi1) * V).ref();
    const volVectorField E2 = (-fvc::grad(phi2) * Vf).ref();
    const volVectorField E3 = (-fvc::grad(phi3) * Vi).ref();
    const volVectorField E4 = (-fvc::grad(phi4) * Vc).ref();
    const interpolationCellPoint<vector> E1Interp(E1);
    const interpolationCellPoint<vector> E2Interp(E2);
    const interpolationCellPoint<vector> E3Interp(E3);
    const interpolationCellPoint<vector> E4Interp(E4);

    // 時間ステップ
    const double timestep = 0.001e-6;
    double time = 0;

    // 結果出力
    const fileName dir = "results";
    if (!isDir(dir)) {
        mkDir(dir);
    }
    if (trajectory) {
        const fileName trajectoryDir = dir / "trajectory";
        if (!isDir(trajectoryDir)) {
            mkDir(trajectoryDir);
        }
    }
    const fileName file = dir / "result";
    OFstream os(file);
    os << "b count fallen r_hit z_hit" << endl;

    Random rng;

    Info << "Starting" << endl;

    for (double b = Bmin.value(); b < Bmax.value() + Bstep.value() / 2; b += Bstep.value()) {
        double bz = Bz.value() * b;
        double br = Br.value() * b;
        if (!corrected) {
            bz = 0;
            br = 0;
        }
        int count = 0;
        int fallen = 0;
        int r_hit = 0;
        int z_hit = 0;
        for (int n = 1; n <= 1000; ++n) {
            scalar angle = 2.0 * M_PI * rng.sample01<scalar>();
            scalar r = 0.0005 + 0.003 * rng.sample01<scalar>();
            scalar theta = rng.sample01<scalar>();
            point p(r * std::cos(angle), -0.008 + r * std::sin(angle), 0.009 + 0.004 * 2.0 * (rng.sample01<scalar>() - 0.5));
            vector v(0, 0, 0);
            List<scalar> posx(0, 0);
            List<scalar> posy(0, 0);
            List<scalar> posz(0, 0);
            // Runge-Kutta法
            for (int i = 0; i < 100000; ++i) {
                // フィラメント
                if (p.x() * p.x() + (p.y() + 0.008) * (p.y() + 0.008) < 0.00025 * 0.00025 &&
                    p.z() > 0.003 && p.z() < 0.013) {
                    fallen++;
                    break;
                }
                // 下のケーブル
                if (p.x() * p.x() + (p.z() - 0.003) * (p.z() - 0.003) < 0.00025 * 0.00025 &&
                    p.y() > -0.003) {
                    fallen++;
                    break;
                }
                // 上のケーブル
                if (p.x() * p.x() + (p.z() - 0.013) * (p.z() - 0.013) < 0.00025 * 0.00025 &&
                    p.y() > -0.004) {
                    fallen++;
                    break;
                }
                // Dee電極
                if ((p.x() + 0.005) * (p.x() + 0.005) + p.y() * p.y() > 0.025 * 0.025 && p.x() < -0.005) {
                    r_hit++;
                    break;
                }
                // 真空容器の壁
                if (p.x() * p.x() + p.y() * p.y() > 0.037 * 0.037) {
                    r_hit++;
                    break;
                }
                // ダミーDee電極
                if (p.x() > 0.002 && p.x() < 0.008 && (p.y() > 0.025 || p.y() < -0.025)) {
                    r_hit++;
                    break;
                }
                // 上下の壁
                if (p.z() < 0.004 || p.z() > 0.013) {
                    z_hit++;
                    break;
                }
                // 検出器
                if (p.z() > 0.004 && p.z() < 0.014 && p.y() > -0.002 && p.y() < 0.002 && p.x() > 0.018 && p.x() < 0.030) {
                    count++;
                    break;
                }
                label cell = mesh.findCell(p);
                if (cell == -1) {
                    cell = mesh.findNearestCell(p);
                }
                double b1 = b - bz * 0.020 + std::abs(p.z() - 0.0115) * bz - std::sqrt(p.x() * p.x() + p.y() * p.y()) * br;
                const vector dv1 = (
                    E1Interp.interpolate(p, cell) * std::sin(2.0 * M_PI * (f.value() * time + theta))
                    + E2Interp.interpolate(p, cell)
                    + E3Interp.interpolate(p, cell)
                    + E4Interp.interpolate(p, cell)
                    + vector(v.y() * b1, -v.x() * b1, 0.0)) * (e / m * timestep).value();

                const point p2 = p + 0.5 * (v + dv1) * timestep;
                double b2 = b - bz * 0.020 + std::abs(p2.z() - 0.0115) * bz - std::sqrt(p2.x() * p2.x() + p2.y() * p2.y()) * br;
                const vector dv2 = (
                    E1Interp.interpolate(p2, cell) * std::sin(2.0 * M_PI * (f.value() * (time + timestep / 2.0) + theta))
                    + E2Interp.interpolate(p2, cell)
                    + E3Interp.interpolate(p2, cell)
                    + E4Interp.interpolate(p2, cell)
                    + vector(v.y() * b2, -v.x() * b2, 0.0)) * (e / m * timestep).value();

                const point p3 = p2 + 0.5 * (v + dv2) * timestep;
                double b3 = b - bz * 0.020 + std::abs(p3.z() - 0.0115) * bz - std::sqrt(p3.x() * p3.x() + p3.y() * p3.y()) * br;
                const vector dv3 = (
                    E1Interp.interpolate(p3, cell) * std::sin(2.0 * M_PI * (f.value() * (time + timestep / 2.0) + theta))
                    + E2Interp.interpolate(p3, cell)
                    + E3Interp.interpolate(p3, cell)
                    + E4Interp.interpolate(p3, cell)
                    + vector(v.y() * b3, -v.x() * b3, 0.0)) * (e / m * timestep).value();

                const point p4 = p3 + (v + dv3) * timestep;
                double b4 = b - bz * 0.020 + std::abs(p4.z() - 0.0115) * bz - std::sqrt(p4.x() * p4.x() + p4.y() * p4.y()) * br;
                const vector dv4 = (
                    E1Interp.interpolate(p4, cell) * std::sin(2.0 * M_PI * (f.value() * (time + timestep) + theta))
                    + E2Interp.interpolate(p4, cell)
                    + E3Interp.interpolate(p4, cell)
                    + E4Interp.interpolate(p4, cell)
                    + vector(v.y() * b4, -v.x() * b4, 0.0)) * (e / m * timestep).value();

                p += (v + (dv1 * 2.0 + dv2 * 2.0 + dv3) / 6.0) * timestep;
                v += (dv1 + dv2 * 2.0 + dv3 * 2.0 + dv4) / 6.0;
                time += timestep;
                posx.append(p.x());
                posy.append(p.y());
                posz.append(p.z());
            }
            if (trajectory) {
                const fileName trajectoryDir = dir / "trajectory" / std::to_string(b);
                if (!isDir(trajectoryDir)) {
                    mkDir(trajectoryDir);
                }
                OFstream os(trajectoryDir / std::to_string(n));
                forAll(posx, i) {
                    os << posx[i] << " " << posy[i] << " " << posz[i] << endl;
                }
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
