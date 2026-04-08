int main() {
    /*for (Real i = 1; i < 30;) {
        test_ic(i);
        i += 0.01;
    }*/

    test_FRS();
    return 0;

    /*double xi[] = {0.001, 0.01, 0.1, 1, 2, 3, 5, 10, 100};
    // double xi[] = {100};
    double sigma[] = {0, 0.01, 0.05, 0.1, 1, 10, 100};

    // double sigma[] = {0};

    for (auto x : xi) {
        for (auto s : sigma) {
            test_reverse_shock(x, s);
        }
    }

    double xi2[] = {0.001, 0.01, 0.1, 1, 5, 10, 100, 1000, 10000};
    Array sigma2 = xt::logspace(std::log10(1e-5), std::log10(100), 100);

    //double sigma2[] = {0, 0, 0, 0, 0, 0};

    for (auto x : xi2) {
        std::ofstream out("rshock-data/crossing-time-" + std::to_string(x) + ".txt");
        for (auto s : sigma2) {
            auto [t_cross, duration, t_dec] = test_reverse_shock(x, s, false);
            out << x << ' ' << s << ' ' << t_cross << ' ' << duration << ' ' << t_dec << std::endl;
        }
    }

    return 0;
}
