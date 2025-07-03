#include <iostream>
#include <cmath>
#include <iomanip>
#include <functional>
#include <vector>
#include <algorithm>
#include <limits>
#include <string>

using namespace std;

// === Wyniki referencyjne z Wolframa ===
vector<double> wolfram_roots1 = { 0.348076736197189 };
vector<double> wolfram_roots2 = { -1.48417, 0.0333518, 1.72233, 3.11798 };
vector<double> wolfram_roots3 = { -2.6176, -2.3569 , -1.6502 , -1.3608 ,-0.963027395181688,
                                -0.709428946764046, -0.282374682169686,
                                -0.0547080458855929, 0.394475324542100, 0.602895578780060,
                                1.069001610280017, 1.262385075266662, 1.7423, 1.9223, 2.4139,
                                2.5850, 3.0841, 3.2469, 3.7559, 3.9100
};

// === Struktury do przechowywania wyników z iteracjami ===
struct MethodResult {
    double root;
    int iterations;
    bool converged;

    MethodResult() : root(NAN), iterations(0), converged(false) {}
    MethodResult(double r, int iter, bool conv) : root(r), iterations(iter), converged(conv) {}
};

// === Struktura dla wiersza tabeli ===
struct TableRow {
    double approx_root;
    string method_name;
    double root;
    double function_value;
    double error;
    int iterations;
    bool converged;

    TableRow(double approx, string name, double r, double fval, double err, int iter, bool conv)
        : approx_root(approx), method_name(name), root(r), function_value(fval), error(err), iterations(iter), converged(conv) {}
};

// === Funkcje ===

double equation1(double x) {
    return log(x + 1) - 1.0 / (x + 3);
}

double equation2(double x) {
    return pow(x, 3) + 30.0 * cos(x) - 1.0 / x;
}

double equation3(double x) {
    const double PI = 3.141592653589793;
    return sin(3 * PI * x) / (x + 2) + 1.0 / (x + 4);
}

// === Pochodne analityczne ===

double derivative1(double x) {
    return 1.0 / (x + 1) + 1.0 / ((x + 3) * (x + 3));
}

double derivative2(double x) {
    return 3 * x * x - 30.0 * sin(x) + 1.0 / (x * x);
}

// === Pochodna numeryczna ===

double numeric_derivative(function<double(double)> f, double x, double h = 1e-6) {
    return (f(x + h) - f(x - h)) / (2 * h);
}

// === Funkcja do znajdowania najbli¿szego wyniku referencyjnego ===

pair<double, double> find_closest_reference(double computed_root, const vector<double>& reference_roots) {
    double min_error = numeric_limits<double>::max();
    double closest_ref = 0;

    for (double ref_root : reference_roots) {
        double error = fabs(computed_root - ref_root);
        if (error < min_error) {
            min_error = error;
            closest_ref = ref_root;
        }
    }

    return make_pair(closest_ref, min_error);
}

// === Funkcje do wyœwietlania tabeli ===

void print_table_header() {
    cout << "+" << string(16, '-') << "+" << string(20, '-') << "+" << string(18, '-')
        << "+" << string(15, '-') << "+" << string(15, '-') << "+" << string(12, '-') << "+" << endl;
    cout << "| " << setw(14) << left << "Przybliz. x"
        << "| " << setw(20) << "Metoda"
        << "| " << setw(16) << "Miejsce zerowe"
        << "| " << setw(13) << "f(x)"
        << "| " << setw(13) << "Blad"
        << "| " << setw(10) << "Iteracje" << "|" << endl;
    cout << "+" << string(16, '-') << "+" << string(20, '-') << "+" << string(18, '-')
        << "+" << string(15, '-') << "+" << string(15, '-') << "+" << string(12, '-') << "+" << endl;
}

void print_table_row(const TableRow& row, bool show_approx = true) {
    if (show_approx) {
        cout << "| " << setw(14) << fixed << setprecision(6) << row.approx_root;
    }
    else {
        cout << "| " << setw(14) << "";
    }

    cout << "| " << setw(20) << left << row.method_name;

    if (row.converged) {
        cout << "| " << setw(16) << fixed << setprecision(8) << row.root
            << "| " << setw(13) << scientific << setprecision(2) << row.function_value
            << "| " << setw(13) << scientific << setprecision(2) << row.error
            << "| " << setw(10) << row.iterations << "|" << endl;
    }
    else {
        cout << "| " << setw(16) << "Nie zbieg³a"
            << "| " << setw(13) << "-"
            << "| " << setw(13) << "-"
            << "| " << setw(10) << "-" << "|" << endl;
    }
}

void print_table_footer() {
    cout << "+" << string(16, '-') << "+" << string(20, '-') << "+" << string(18, '-')
        << "+" << string(15, '-') << "+" << string(15, '-') << "+" << string(12, '-') << "+" << endl;
}

void print_reference_roots(const vector<double>& roots, const string& equation_name) {
    cout << "\n=== " << equation_name << " ===" << endl;
    cout << "Wyniki referencyjne z Wolframa: ";
    for (size_t i = 0; i < roots.size(); ++i) {
        cout << fixed << setprecision(6) << roots[i];
        if (i < roots.size() - 1) cout << ", ";
    }
    cout << endl << endl;
}

// === Metoda bisekcji ===

MethodResult bisection(function<double(double)> f, double a, double b, double epsilon = 1e-8, int maxIterations = 1000) {
    if (f(a) * f(b) >= 0) {
        return MethodResult();
    }

    double c;
    int iteration = 0;
    while ((b - a) / 2.0 > epsilon && iteration < maxIterations) {
        iteration++;
        c = (a + b) / 2.0;
        double fc = f(c);

        if (fabs(fc) < epsilon)
            return MethodResult(c, iteration, true);

        if (f(a) * fc < 0)
            b = c;
        else
            a = c;
    }

    return MethodResult((a + b) / 2.0, iteration, true);
}

// === Metoda Newtona (pochodna numeryczna) ===

MethodResult newton_numeric(function<double(double)> f, double x0, double epsilon = 1e-8, int maxIterations = 1000) {
    double x = x0;
    for (int i = 0; i < maxIterations; ++i) {
        double fx = f(x);
        double dfx = numeric_derivative(f, x);

        if (fabs(dfx) < 1e-14) {
            return MethodResult();
        }

        double x_new = x - fx / dfx;
        if (fabs(x_new - x) < epsilon)
            return MethodResult(x_new, i + 1, true);

        x = x_new;
    }
    return MethodResult();
}

// === Metoda Newtona (pochodna analityczna) ===

MethodResult newton_analytic(function<double(double)> f, function<double(double)> df, double x0, double epsilon = 1e-8, int maxIterations = 1000) {
    double x = x0;
    for (int i = 0; i < maxIterations; ++i) {
        double fx = f(x);
        double dfx = df(x);

        if (fabs(dfx) < 1e-14) {
            return MethodResult();
        }

        double x_new = x - fx / dfx;
        if (fabs(x_new - x) < epsilon)
            return MethodResult(x_new, i + 1, true);

        x = x_new;
    }
    return MethodResult();
}

// === Metoda siecznych ===

MethodResult secant_method(function<double(double)> f, double x0, double x1, double epsilon = 1e-8, int maxIterations = 1000) {
    double f0 = f(x0);
    double f1 = f(x1);

    for (int i = 0; i < maxIterations; ++i) {
        if (fabs(f1 - f0) < 1e-14) {
            return MethodResult();
        }

        double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
        if (fabs(x2 - x1) < epsilon)
            return MethodResult(x2, i + 1, true);

        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f(x1);
    }

    return MethodResult();
}

// === Funkcja do szukania wszystkich pierwiastków metod¹ bisekcji w zadanym przedziale ===

vector<double> find_all_roots(function<double(double)> f, double start, double end, double step = 0.05, double epsilon = 1e-8) {
    vector<double> roots;
    double a = start;
    double b = a + step;

    while (b <= end) {
        if (isnan(f(a)) || isnan(f(b))) {
            a = b;
            b = a + step;
            continue;
        }

        double fa = f(a);
        double fb = f(b);

        if (fa * fb < 0) {
            MethodResult result = bisection(f, a, b, epsilon);
            if (result.converged) {
                bool is_new = true;
                for (double r : roots) {
                    if (fabs(r - result.root) < epsilon * 10) {
                        is_new = false;
                        break;
                    }
                }
                if (is_new) {
                    roots.push_back(result.root);
                }
            }
        }

        a = b;
        b = a + step;
    }

    return roots;
}

// === Funkcja do analizy jednego pierwiastka ===

vector<TableRow> analyze_root(function<double(double)> f, function<double(double)> df, double root,
    const vector<double>& reference_roots, bool has_analytic_derivative = true) {
    vector<TableRow> results;

    // Bisekcja
    MethodResult bis_result = bisection(f, root - 0.1, root + 0.1);
    if (bis_result.converged) {
        auto bis_error = find_closest_reference(bis_result.root, reference_roots);
        results.emplace_back(root, "Bisekcja", bis_result.root, f(bis_result.root),
            bis_error.second, bis_result.iterations, true);
    }
    else {
        results.emplace_back(root, "Bisekcja", 0, 0, 0, 0, false);
    }

    // Newton analityczna (jeœli dostêpna)
    if (has_analytic_derivative) {
        MethodResult newt_a_result = newton_analytic(f, df, root);
        if (newt_a_result.converged) {
            auto newt_a_error = find_closest_reference(newt_a_result.root, reference_roots);
            results.emplace_back(root, "Newton (analityczna)", newt_a_result.root, f(newt_a_result.root),
                newt_a_error.second, newt_a_result.iterations, true);
        }
        else {
            results.emplace_back(root, "Newton (analityczna)", 0, 0, 0, 0, false);
        }
    }

    // Newton numeryczna
    MethodResult newt_n_result = newton_numeric(f, root);
    if (newt_n_result.converged) {
        auto newt_n_error = find_closest_reference(newt_n_result.root, reference_roots);
        results.emplace_back(root, "Newton (numeryczna)", newt_n_result.root, f(newt_n_result.root),
            newt_n_error.second, newt_n_result.iterations, true);
    }
    else {
        results.emplace_back(root, "Newton (numeryczna)", 0, 0, 0, 0, false);
    }

    // Sieczne
    MethodResult sec_result = secant_method(f, root - 0.1, root + 0.1);
    if (sec_result.converged) {
        auto sec_error = find_closest_reference(sec_result.root, reference_roots);
        results.emplace_back(root, "Sieczne", sec_result.root, f(sec_result.root),
            sec_error.second, sec_result.iterations, true);
    }
    else {
        results.emplace_back(root, "Sieczne", 0, 0, 0, 0, false);
    }

    return results;
}

// === MAIN ===

int main() {
    double start = -3.0;
    double end = 4.0;

    // --- Funkcja 1 ---
    print_reference_roots(wolfram_roots1, "Funkcja 1: ln(x+1) - 1/(x+3)");

    auto roots1 = find_all_roots(equation1, start, end);
    if (roots1.empty()) {
        cout << "Nie znaleziono miejsc zerowych." << endl << endl;
    }
    else {
        vector<TableRow> all_results1;
        for (auto root : roots1) {
            auto root_results = analyze_root(equation1, derivative1, root, wolfram_roots1, true);
            all_results1.insert(all_results1.end(), root_results.begin(), root_results.end());
        }

        print_table_header();
        double last_approx = -999;
        for (const auto& row : all_results1) {
            bool show_approx = (row.approx_root != last_approx);
            print_table_row(row, show_approx);
            last_approx = row.approx_root;
        }
        print_table_footer();
        cout << endl;
    }

    // --- Funkcja 2 ---
    print_reference_roots(wolfram_roots2, "Funkcja 2: x^3 + 30*cos(x) - 1/x");

    auto roots2 = find_all_roots(equation2, start, end);
    if (roots2.empty()) {
        cout << "Nie znaleziono miejsc zerowych." << endl << endl;
    }
    else {
        vector<TableRow> all_results2;
        for (auto root : roots2) {
            auto root_results = analyze_root(equation2, derivative2, root, wolfram_roots2, true);
            all_results2.insert(all_results2.end(), root_results.begin(), root_results.end());
        }

        print_table_header();
        double last_approx = -999;
        for (const auto& row : all_results2) {
            bool show_approx = (row.approx_root != last_approx);
            print_table_row(row, show_approx);
            last_approx = row.approx_root;
        }
        print_table_footer();
        cout << endl;
    }

    // --- Funkcja 3 ---
    print_reference_roots(wolfram_roots3, "Funkcja 3: sin(3PIx)/(x+2) + 1/(x+4)");

    auto roots3 = find_all_roots(equation3, start, end);
    if (roots3.empty()) {
        cout << "Nie znaleziono miejsc zerowych." << endl << endl;
    }
    else {
        vector<TableRow> all_results3;
        for (auto root : roots3) {
            auto root_results = analyze_root(equation3, nullptr, root, wolfram_roots3, false);
            all_results3.insert(all_results3.end(), root_results.begin(), root_results.end());
        }

        print_table_header();
        double last_approx = -999;
        for (const auto& row : all_results3) {
            bool show_approx = (row.approx_root != last_approx);
            print_table_row(row, show_approx);
            last_approx = row.approx_root;
        }
        print_table_footer();
        cout << endl;
    }

    return 0;
}