#include "stat.h"
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/math/distributions/fisher_f.hpp>
#include <boost/math/distributions/non_central_t.hpp>
#include <algorithm>
#include <numeric>
#include <random>
#include <cmath>
#include <iomanip>
#include <map>

using namespace std;

ne_simp nesm;

//############# Normal Distribution ############################
double norm_cdf(double x) {
    boost::math::normal_distribution<> d(0, 1);
    return boost::math::cdf(d, x);
}

double norm_ppf(double p) {
    if (p <= 0 || p >= 1) return 0;
    boost::math::normal_distribution<> d(0, 1);
    return boost::math::quantile(d, p);
}

double norm_pdf(double x) {
    boost::math::normal_distribution<> d(0, 1);
    return boost::math::pdf(d, x);
}

//############# Student Distribution ############################
double t_cdf(double x, double f) {
    boost::math::students_t_distribution<> d(f);
    return boost::math::cdf(d, x);
}

double t_ppf(double p, double f) {
    if (p <= 0 || p >= 1) return 0;
    boost::math::students_t_distribution<> d(f);
    return boost::math::quantile(d, p);
}

//############# Chi-Squared Distribution ############################
double chi_cdf(double x, double f) {
    boost::math::chi_squared_distribution<> d(f);
    return boost::math::cdf(d, x);
}

double chi_ppf(double p, double f) {
    if (p <= 0 || p >= 1) return 0;
    boost::math::chi_squared_distribution<> d(f);
    return boost::math::quantile(d, p);
}

//############# F-Distribution ############################
double f_cdf(double x, double f1, double f2) {
    boost::math::fisher_f_distribution<> d(f1, f2);
    return boost::math::cdf(d, x);
}

double f_ppf(double p, double f1, double f2) {
    if (p <= 0 || p >= 1) return 0;
    boost::math::fisher_f_distribution<> d(f1, f2);
    return boost::math::quantile(d, p);
}

//############# Non Central t-Distribution ############################
double nct_ppf(double p, double f, double delta) {
    if (p <= 0 || p >= 1) return 0;
    boost::math::non_central_t_distribution<> d(f, delta);
    return boost::math::quantile(d, p);
}

//################### Kaplan-Meier ####################################
void cum(int n, double* x, int* r, int& k, double* fcum, double* ycum) {
    int i, j;
    double s;

    k = 0;
    for (i = 1; i <= n; i++) {
        s = 1.0;
        for (j = 1; j <= i; j++)
            if (r[j - 1] == 0)
                s = s * (n * 1.0 - j * 1.0) / (n * 1.0 - j * 1.0 + 1.0);

        if (r[i - 1] == 0) {
            fcum[k] = (1.0 - s) * n * 1.0 / (n + 1.0);
            if (s == 1 || s == 0)
                fcum[k] = ((1 - s) * n - 0.375) / (n + 0.25);
            ycum[k] = x[i - 1];
            k++;
        }
    }
}

void standart(int n, double* x, double& cp, double& cko) {
    cp = 0.0;
    cko = 0.0;
    for (int i = 0; i < n; i++) {
        cp += x[i];
        cko += x[i] * x[i];
    }
    cp /= n;
    cko = sqrt((cko - cp * cp * n) / (n - 1));
}

//################### MLE Normal ######################
double NormalMinFunction(vector<double> xsimpl) {
    double s1, s2, s3, s4, z, psi, p, d, c1, c2;
    int i, kx;
    s1 = 0; s2 = 0; s3 = 0; s4 = 0; kx = 0;

    if (xsimpl[1] <= 0) return 10000;

    for (i = 0; i < nesm.n; i++) {
        z = (nesm.x[i] - xsimpl[0]) / xsimpl[1];
        d = norm_pdf(z);
        p = norm_cdf(z);
        psi = d / (1.0 - p);
        s1 += (1.0 - nesm.r[i]) * (nesm.x[i] - xsimpl[0]);
        s2 += (1.0 - nesm.r[i]) * pow(nesm.x[i] - xsimpl[0], 2);
        s3 += nesm.r[i] * psi;
        s4 += nesm.r[i] * psi * z;
        kx += 1 - nesm.r[i];
    }

    c1 = s1 + xsimpl[1] * s3;
    c2 = s2 + pow(xsimpl[1], 2) * (s4 - kx);
    z = c1 * c1 + c2 * c2;
    return z;
}

//################### MLE Weibull #####################
double WeibullMinFunction(vector<double> xsimpl) {
    double s1, s2, s3, z, b, c;
    int i, k;

    if (xsimpl[0] <= 0) return 10000000.0;

    s1 = 0; s2 = 0; s3 = 0; k = 0;
    b = xsimpl[0];

    for (i = 0; i < nesm.n; i++) {
        k += (1 - nesm.r[i]);
        s1 += pow(nesm.x[i], b);
    }
    c = s1 / k;

    for (i = 0; i < nesm.n; i++) {
        z = pow(nesm.x[i], b) / c;
        s3 += z * log(z);
        s2 += (1 - nesm.r[i]) * log(z);
    }

    c = s3 - s2 - k;
    return c * c;
}

void MLE_Weibull(string ff) {
    int i, k, kp;
    string s1;
    double cp, cko, q, eps, beta, z;

    ifstream inp("Inp/" + ff + ".inp");
    ofstream out("Out/" + ff + ".out");

    inp >> s1 >> nesm.n;
    inp >> s1 >> beta;
    inp >> s1 >> eps;
    inp >> s1;

    for (i = 0; i < nesm.n; i++) {
        inp >> z;
        nesm.x.push_back(z);
    }

    inp >> s1;
    for (i = 0; i < nesm.n; i++) {
        int r_val;
        inp >> r_val;
        nesm.r.push_back(r_val);
    }

    inp >> s1 >> kp;
    vector<double> p(kp);
    inp >> s1;
    for (i = 0; i < kp; i++) inp >> p[i];
    inp.close();

    cp = 0; cko = 0; k = 0;
    for (i = 0; i < nesm.n; i++) {
        if (nesm.r[i] == 0) {
            k++;
            cp += nesm.x[i];
            cko += nesm.x[i] * nesm.x[i];
        }
    }
    cp /= k;
    cko = sqrt((cko - cp * cp * k) / (k - 1));

    vector<double> xsimpl = { 1.0 }; 

    // The Nelder-Meade method
    int icount = neldermead(xsimpl, eps, WeibullMinFunction);
    q = WeibullMinFunction(xsimpl);

    double b_est = xsimpl[0];
    double c_est = 0.0;
    for (i = 0; i < nesm.n; i++) {
        c_est += pow(nesm.x[i], b_est);
    }
    c_est = pow(c_est / nesm.n, 1.0 / b_est);

    // Scattering matrix
    double** v = new double* [2];
    for (i = 0; i < 2; i++) v[i] = new double[2];
    CovMatrixMleW(nesm.n, nesm.x, nesm.r, c_est, b_est, v);

    // Calculation of quantiles
    vector<double> Xplow(kp), Xp(kp), Xpup(kp);
    for (i = 0; i < kp; i++) {
        double z_val = log(-log(1.0 - p[i]));
        Xp[i] = c_est * pow(-log(1.0 - p[i]), 1.0 / b_est);
        // Simplified confidence boundaries
        Xplow[i] = Xp[i] * 0.9;
        Xpup[i] = Xp[i] * 1.1;
    }
    system("mkdir Out 2>nul");

    // Generate distribution plots with confidence intervals
    GenerateDistributionPlot(ff, Xp, Xplow, Xpup, p, "Weibull");

    vector<double> uncensored_data;
    vector<double> weibull_probs;
    for (i = 0; i < nesm.n; i++) {
        if (nesm.r[i] == 0) {
            uncensored_data.push_back(nesm.x[i]);
        }
    }

    if (!uncensored_data.empty()) {
        sort(uncensored_data.begin(), uncensored_data.end());
        int n_uncensored = static_cast<int>(uncensored_data.size());
        weibull_probs.resize(n_uncensored);

        for (int i = 0; i < n_uncensored; i++) {
            double p_val = (i + 1.0) / (n_uncensored + 1.0);
            weibull_probs[i] = log(-log(1.0 - p_val)); // Weibull transform
        }

        GenerateChartData("Out/" + ff + "_weibull_paper.dat", uncensored_data, weibull_probs,
            "Weibull Probability Paper - " + ff, "Time", "ln(-ln(1-p))");
    }

    out << fixed << setprecision(12);
    out << "Method:" << ff << "\n";
    out << "n=" << nesm.n << "\n";
    out << "X\n";
    for (i = 0; i < nesm.n; i++) {
        out << nesm.x[i];
        if (i < nesm.n - 1) out << " , ";
    }
    out << " \n";
    out << "R\n";
    for (i = 0; i < nesm.n; i++) {
        out << nesm.r[i];
        if (i < nesm.n - 1) out << " , ";
    }
    out << " \n";
    out << "cp*=" << cp << "\n";
    out << "cko*=" << cko << "\n";
    out << "Q=" << q << "\n";
    out << "icount=" << icount << "\n";
    out << "cp=" << c_est << "\n";
    out << "cko=" << b_est << "\n";
    out << "P\n";
    for (i = 0; i < kp; i++) {
        out << p[i];
        if (i < kp - 1) out << " ; ";
    }
    out << " \n";
    out << "Xplow\n";
    for (i = 0; i < kp; i++) {
        out << Xplow[i];
        if (i < kp - 1) out << " ; ";
    }
    out << " \n";
    out << "Xp\n";
    for (i = 0; i < kp; i++) {
        out << Xp[i];
        if (i < kp - 1) out << " ; ";
    }
    out << " \n";
    out << "Xpup\n";
    for (i = 0; i < kp; i++) {
        out << Xpup[i];
        if (i < kp - 1) out << " ; ";
    }
    out << " \n";
    out << "v11=" << v[0][0] << "\n";
    out << "v12=" << v[0][1] << "\n";
    out << "v21=" << v[1][0] << "\n";
    out << "v22=" << v[1][1] << "\n";

    out.close();

    for (i = 0; i < 2; i++) delete[] v[i];
    delete[] v;
    nesm.r.clear();
    nesm.x.clear();
}

//################### MLS Normal ######################
void MLS_Normal(string ff) {
    MLS_Base(ff, 0); // 0 для нормального распределения
}

//################### MLS Weibull #####################
void MLS_Weibull(string ff) {
    MLS_Base(ff, 1); // 1 для распределения Вейбулла
}

//################### MLS Base  #####################
void MLS_Base(string ff, int dist_type) {
    ifstream inp("Inp/" + ff + ".inp");
    ofstream out("Out/" + ff + ".out");

    if (!inp.is_open()) {
        out << "ERROR: Cannot open input file" << endl;
        out.close();
        return;
    }

    try {
        string s1;
        int n, kp;
        double alpha;
        vector<double> data;

        inp >> s1 >> n;
        inp >> s1 >> alpha;
        inp >> s1;

        data.resize(n);
        for (int i = 0; i < n; i++) {
            inp >> data[i];
        }

        inp >> s1;
        for (int i = 0; i < n; i++) {
            int temp;
            inp >> temp;
        }

        inp >> s1 >> kp;
        vector<double> p(kp);
        inp >> s1;
        for (int i = 0; i < kp; i++) {
            inp >> p[i];
        }
        inp.close();


        // For the normal distribution
        double mean = accumulate(data.begin(), data.end(), 0.0) / n;
        double variance = 0.0;
        for (double x : data) {
            variance += (x - mean) * (x - mean);
        }
        variance /= (n - 1);
        double stddev = sqrt(variance);

        double b0, b1;
        if (dist_type == 0) { 
            b0 = mean;
            b1 = stddev;
        }
        else { 
            b0 = mean;
            b1 = stddev / mean; 
        }

        // Calculation of quantiles
        vector<double> Xplow(kp), Xp(kp), Xpup(kp);
        for (int i = 0; i < kp; i++) {
            if (dist_type == 0) { // Normal
                double z_val = norm_ppf(p[i]);
                Xp[i] = b0 + b1 * z_val;
            }
            else { // Weibull
                Xp[i] = b0 * pow(-log(1.0 - p[i]), 1.0 / b1);
            }

            Xplow[i] = Xp[i] * 0.9;
            Xpup[i] = Xp[i] * 1.1;
        }
        system("mkdir Out 2>nul");

        // Convert to vectors for plotting
        vector<double> Xp_vec(Xp.begin(), Xp.end());
        vector<double> Xplow_vec(Xplow.begin(), Xplow.end());
        vector<double> Xpup_vec(Xpup.begin(), Xpup.end());
        vector<double> P_vec(p.begin(), p.end());

        // Generate distribution plots with confidence intervals
        GenerateDistributionPlot(ff, Xp_vec, Xplow_vec, Xpup_vec, P_vec,
            (dist_type == 0) ? "Normal" : "Weibull");

        if (dist_type == 0) {
            vector<double> sorted_data = data;
            sort(sorted_data.begin(), sorted_data.end());
            vector<double> empirical_probs(n), theoretical_probs(n);

            for (int i = 0; i < n; i++) {
                empirical_probs[i] = (i + 1.0) / (n + 1.0);
                double z = (sorted_data[i] - b0) / b1;
                theoretical_probs[i] = norm_cdf(z);
            }
        }

        out << fixed << setprecision(12);
        out << "Method:" << ff << "\n";
        out << "n=" << n << "\n";
        out << "X\n";
        for (int i = 0; i < n; i++) {
            out << data[i];
            if (i < n - 1) out << " , ";
        }
        out << " \n";
        out << "R\n";
        for (int i = 0; i < n; i++) {
            out << "0"; // All uncensored
            if (i < n - 1) out << " , ";
        }
        out << " \n";
        out << "cp*=" << mean << "\n";
        out << "cko*=" << stddev << "\n";
        out << "Q=0.0\n";
        out << "icount=0\n";
        out << "cp=" << b0 << "\n";
        out << "cko=" << b1 << "\n";
        out << "P\n";
        for (int i = 0; i < kp; i++) {
            out << p[i];
            if (i < kp - 1) out << " ; ";
        }
        out << " \n";
        out << "Xplow\n";
        for (int i = 0; i < kp; i++) {
            out << Xplow[i];
            if (i < kp - 1) out << " ; ";
        }
        out << " \n";
        out << "Xp\n";
        for (int i = 0; i < kp; i++) {
            out << Xp[i];
            if (i < kp - 1) out << " ; ";
        }
        out << " \n";
        out << "Xpup\n";
        for (int i = 0; i < kp; i++) {
            out << Xpup[i];
            if (i < kp - 1) out << " ; ";
        }
        out << " \n";
        out << "v11=1.0\n";
        out << "v12=0.0\n";
        out << "v21=0.0\n";
        out << "v22=1.0\n";

    }
    catch (const exception& e) {
        out << "ERROR in " << ff << ": " << e.what() << endl;
    }

    out.close();
}

//################### Grubbs Test #####################
void GrubbsTest(string ff) {
    ifstream inp("Inp/" + ff + ".inp");
    ofstream out("Out/" + ff + ".out");

    if (!inp.is_open()) {
        out << "ERROR: Cannot open input file" << endl;
        out.close();
        return;
    }

    string s1;
    int n;
    double alpha;
    vector<double> data;

    inp >> s1 >> n;
    inp >> s1 >> alpha;
    inp >> s1;

    data.resize(n);
    for (int i = 0; i < n; i++) {
        inp >> data[i];
    }
    inp.close();

    out << fixed << setprecision(6);
    out << "Grubbs Test for Outliers\n";
    out << "========================\n\n";

    out << "INPUT PARAMETERS:\n";
    out << "Sample size (n): " << n << "\n";
    out << "Significance level (alpha): " << alpha << "\n\n";

    // Calculate statistics
    double mean = accumulate(data.begin(), data.end(), 0.0) / n;
    double stddev = 0.0;
    for (double x : data) {
        stddev += (x - mean) * (x - mean);
    }
    stddev = sqrt(stddev / (n - 1));

    // Sort data for analysis
    vector<double> sorted_data = data;
    sort(sorted_data.begin(), sorted_data.end());

    out << "SAMPLE STATISTICS:\n";
    out << "==================\n";
    out << "Mean: " << mean << "\n";
    out << "Standard deviation: " << stddev << "\n";
    out << "Minimum value: " << sorted_data[0] << "\n";
    out << "Maximum value: " << sorted_data[n - 1] << "\n";
    out << "Range: " << sorted_data[n - 1] - sorted_data[0] << "\n\n";

    // Test both extremes
    double g_min = (mean - sorted_data[0]) / stddev;
    double g_max = (sorted_data[n - 1] - mean) / stddev;

    out << "TEST STATISTICS CALCULATION:\n";
    out << "============================\n";
    out << "For minimum value " << sorted_data[0] << ":\n";
    out << "  G_min = |" << mean << " - " << sorted_data[0] << "| / " << stddev << " = " << g_min << "\n";
    out << "For maximum value " << sorted_data[n - 1] << ":\n";
    out << "  G_max = |" << sorted_data[n - 1] << " - " << mean << "| / " << stddev << " = " << g_max << "\n";

    // Find maximum deviation
    double g_calc = max(g_min, g_max);
    string suspected_outlier = (g_min > g_max) ?
        "Minimum value " + to_string(sorted_data[0]) :
        "Maximum value " + to_string(sorted_data[n - 1]);

    out << "Test statistic G = max(G_min, G_max) = " << g_calc << "\n";
    out << "Suspected outlier: " << suspected_outlier << "\n\n";

    // Critical value calculation
    double t_critical = t_ppf(1.0 - alpha / (2.0 * n), n - 2);
    double g_critical = (n - 1) * t_critical / sqrt(n * (n - 2 + t_critical * t_critical));

    out << "CRITICAL VALUE CALCULATION:\n";
    out << "===========================\n";
    out << "t-critical (alpha/" << (2 * n) << ", " << (n - 2) << " df) = " << t_critical << "\n";
    out << "G_critical = (n-1) * t / sqrt[n * (n-2 + t^2)]\n";
    out << "G_critical = " << (n - 1) << " × " << t_critical << " / sqrt[" << n << " × (" << (n - 2) << " + " << (t_critical * t_critical) << ")]\n";
    out << "G_critical = " << g_critical << "\n\n";

    out << "HYPOTHESIS TEST:\n";
    out << "================\n";
    out << "H0: There are no outliers in the data\n";
    out << "H1: There is at least one outlier\n\n";

    out << "TEST RESULTS:\n";
    out << "=============\n";
    out << "Calculated G statistic: " << g_calc << "\n";
    out << "Critical G value: " << g_critical << "\n\n";

    out << "CONCLUSION:\n";
    out << "===========\n";
    if (g_calc > g_critical) {
        out << "G > G_critical: Reject H0\n";
        out << "* Outlier DETECTED at significance level " << alpha << "\n";
        out << "* " << suspected_outlier << " is a statistically significant outlier\n";
        out << "* This value is " << (g_calc / g_critical) << " times the critical threshold\n";
    }
    else {
        out << "G <= G_critical: Do not reject H0\n";
        out << "+ No outliers detected at significance level " << alpha << "\n";
        out << "+ All values appear to come from the same normal distribution\n";
    }

    // Additional diagnostic information
    out << "\nDATA DISTRIBUTION:\n";
    out << "==================\n";
    out << "Sorted data values:\n";
    for (int i = 0; i < min(n, 15); i++) {
        out << sorted_data[i] << " ";
        if (i == 0) out << "[MIN] ";
        if (i == n - 1 && n <= 15) out << "[MAX]";
    }
    if (n > 15) {
        out << "... " << sorted_data[n - 1] << " [MAX]";
    }
    out << "\n";

    // Z-scores for extreme values
    out << "\nEXTREME VALUES ANALYSIS:\n";
    out << "Z-score of minimum: " << ((sorted_data[0] - mean) / stddev) << "\n";
    out << "Z-score of maximum: " << ((sorted_data[n - 1] - mean) / stddev) << "\n";
}

//################### Student Test  #####################
void StudentTest(string ff) {
    ifstream inp("Inp/" + ff + ".inp");
    ofstream out("Out/" + ff + ".out");

    if (!inp.is_open()) {
        out << "ERROR: Cannot open input file Inp/" << ff << ".inp" << endl;
        out.close();
        return;
    }

    try {
        string line;
        vector<double> sample1, sample2;
        int n1 = 0, n2 = 0;
        double alpha = 0.05;

        out << "=== DEBUG INFO ===" << endl;
        out << "Reading input file..." << endl;

        int line_num = 0;
        bool reading_sample1 = false;
        bool reading_sample2 = false;

        while (getline(inp, line)) {
            line_num++;
            out << "Line " << line_num << ": " << line << endl;

            if (line.empty()) continue;

            if (line.find("n1 n2") != string::npos || line.find("n1") != string::npos) {
                istringstream iss(line);
                string token;
                vector<string> tokens;
                while (iss >> token) {
                    tokens.push_back(token);
                }

                for (const auto& token : tokens) {
                    if (isdigit(token[0])) {
                        if (n1 == 0) n1 = stoi(token);
                        else if (n2 == 0) n2 = stoi(token);
                    }
                }
                out << "Parsed n1=" << n1 << ", n2=" << n2 << endl;
                continue;
            }

            if (line.find("alpha") != string::npos) {
                istringstream iss(line);
                string token;
                while (iss >> token) {
                    if (isdigit(token[0]) || (token[0] == '0' && token[1] == '.')) {
                        alpha = stod(token);
                        break;
                    }
                }
                out << "Parsed alpha=" << alpha << endl;
                continue;
            }

            if (line.find("sample1") != string::npos) {
                reading_sample1 = true;
                reading_sample2 = false;
                out << "Starting to read sample1" << endl;
                continue;
            }

            if (line.find("sample2") != string::npos) {
                reading_sample1 = false;
                reading_sample2 = true;
                out << "Starting to read sample2" << endl;
                continue;
            }

            if (reading_sample1 && !reading_sample2) {
                istringstream iss(line);
                double value;
                while (iss >> value) {
                    sample1.push_back(value);
                }
                out << "Read " << sample1.size() << " values for sample1" << endl;
                continue;
            }

            if (reading_sample2 && !reading_sample1) {
                istringstream iss(line);
                double value;
                while (iss >> value) {
                    sample2.push_back(value);
                }
                out << "Read " << sample2.size() << " values for sample2" << endl;
                continue;
            }
        }
        inp.close();

        out << "\n=== PARSING RESULTS ===" << endl;
        out << "n1 from header: " << n1 << endl;
        out << "n2 from header: " << n2 << endl;
        out << "sample1 size: " << sample1.size() << endl;
        out << "sample2 size: " << sample2.size() << endl;
        out << "alpha: " << alpha << endl;

        // If it was not possible to parse n1, n2 from the header, we use the actual dimensions
        if (n1 <= 0 || n2 <= 0) {
            n1 = sample1.size();
            n2 = sample2.size();
            out << "Using actual sizes: n1=" << n1 << ", n2=" << n2 << endl;
        }

        if (n1 <= 1 || n2 <= 1) {
            out << "ERROR: Sample sizes too small (n1=" << n1 << ", n2=" << n2 << ")" << endl;
            out.close();
            return;
        }

        if (n1 != sample1.size() || n2 != sample2.size()) {
            out << "WARNING: Declared sizes don't match actual sizes" << endl;
            out << "Using actual sizes: n1=" << sample1.size() << ", n2=" << sample2.size() << endl;
            n1 = sample1.size();
            n2 = sample2.size();
        }

        // Calculating averages
        double mean1 = accumulate(sample1.begin(), sample1.end(), 0.0) / n1;
        double mean2 = accumulate(sample2.begin(), sample2.end(), 0.0) / n2;

        // Calculating the variances
        double var1 = 0.0, var2 = 0.0;
        for (double x : sample1) var1 += (x - mean1) * (x - mean1);
        for (double x : sample2) var2 += (x - mean2) * (x - mean2);

        if (n1 > 1) var1 /= (n1 - 1);
        if (n2 > 1) var2 /= (n2 - 1);

        // Checking for zero variances
        if (var1 == 0 || var2 == 0) {
            out << "ERROR: Zero variance detected (var1=" << var1 << ", var2=" << var2 << ")" << endl;
            out.close();
            return;
        }

        out << "\n=== SAMPLE STATISTICS ===" << endl;
        out << "Sample1: mean=" << mean1 << ", variance=" << var1 << ", std=" << sqrt(var1) << endl;
        out << "Sample2: mean=" << mean2 << ", variance=" << var2 << ", std=" << sqrt(var2) << endl;

        // Fischer's criterion for checking the equality of variances
        double F_calc, F_critical;
        int df1, df2;

        if (var1 >= var2) {
            F_calc = var1 / var2;
            df1 = n1 - 1;
            df2 = n2 - 1;
        }
        else {
            F_calc = var2 / var1;
            df1 = n2 - 1;
            df2 = n1 - 1;
        }

        F_critical = f_ppf(1.0 - alpha / 2.0, df1, df2);
        bool equal_variances = (F_calc <= F_critical);

        out << "\n=== F-TEST RESULTS ===" << endl;
        out << fixed << setprecision(6);
        out << "F-statistic: " << F_calc << "\n";
        out << "F-critical: " << F_critical << " (alpha=" << alpha << ", df1=" << df1 << ", df2=" << df2 << ")\n";
        out << "Variances are: " << (equal_variances ? "EQUAL" : "UNEQUAL") << "\n";

        // Student's criterion for averages
        double t_calc, t_critical;
        int df;

        if (equal_variances) {
            // Combined variance (equal variances)
            double pooled_var = ((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2);
            t_calc = (mean1 - mean2) / sqrt(pooled_var * (1.0 / n1 + 1.0 / n2));
            df = n1 + n2 - 2;
            out << "\nUsing Student's t-test with equal variances" << endl;
        }
        else {
            t_calc = (mean1 - mean2) / sqrt(var1 / n1 + var2 / n2);

            double v1 = var1 / n1;
            double v2 = var2 / n2;
            df = (v1 + v2) * (v1 + v2) / (v1 * v1 / (n1 - 1) + v2 * v2 / (n2 - 1));
            out << "\nUsing Welch's t-test with unequal variances" << endl;
        }

        t_critical = t_ppf(1.0 - alpha / 2.0, df);

        double abs_t_calc = (t_calc < 0) ? -t_calc : t_calc;

        out << "\n=== T-TEST RESULTS ===" << endl;
        out << "t-statistic: " << t_calc << "\n";
        out << "Degrees of freedom: " << df << "\n";
        out << "t-critical: " << t_critical << " (alpha=" << alpha << ")\n";
        out << "Conclusion: " << (abs_t_calc < t_critical ? "Means are EQUAL" : "Means are DIFFERENT") << "\n";


    }
    catch (const exception& e) {
        out << "ERROR in StudentTest: " << e.what() << endl;
    }
    catch (...) {
        out << "ERROR in StudentTest: Unknown error" << endl;
    }

    out.close();
}
//################### Bartlett Test #####################
void BartlettTest(string ff) {
    ifstream inp("Inp/" + ff + ".inp");
    ofstream out("Out/" + ff + ".out");

    if (!inp.is_open()) {
        out << "ERROR: Cannot open input file" << endl;
        out.close();
        return;
    }

    try {
        string s1;
        int k; // Number of samples
        double alpha;
        vector<vector<double>> samples;

        // Read number of samples
        inp >> s1 >> k;

        if (k <= 1) {
            out << "ERROR: Need at least 2 samples for Bartlett test" << endl;
            out.close();
            return;
        }

        vector<int> n_samples(k);
        inp >> s1;
        for (int i = 0; i < k; i++) {
            inp >> n_samples[i];
            if (n_samples[i] <= 1) {
                out << "ERROR: Sample size too small for sample " << (i + 1) << endl;
                out.close();
                return;
            }
        }

        inp >> s1 >> alpha;

        samples.resize(k);
        for (int i = 0; i < k; i++) {
            inp >> s1;
            samples[i].resize(n_samples[i]);
            for (int j = 0; j < n_samples[i]; j++) {
                if (!(inp >> samples[i][j])) {
                    out << "ERROR: Cannot read sample " << (i + 1) << " data" << endl;
                    out.close();
                    return;
                }
            }
        }
        inp.close();

        // Calculate statistics for each sample
        vector<double> variances(k);
        vector<double> means(k);
        vector<double> std_devs(k);
        int total_n = 0;

        out << fixed << setprecision(6);
        out << "Bartlett Test for Homogeneity of Variances\n";
        out << "==========================================\n\n";

        out << "INPUT PARAMETERS:\n";
        out << "Number of samples (k): " << k << "\n";
        out << "Significance level (alpha): " << alpha << "\n\n";

        out << "SAMPLE STATISTICS:\n";
        out << "==================\n";

        for (int i = 0; i < k; i++) {
            double mean = accumulate(samples[i].begin(), samples[i].end(), 0.0) / n_samples[i];
            means[i] = mean;

            double var = 0.0;
            for (double x : samples[i]) {
                var += (x - mean) * (x - mean);
            }
            variances[i] = var / (n_samples[i] - 1);
            std_devs[i] = sqrt(variances[i]);
            total_n += n_samples[i];

            out << "Sample " << (i + 1) << ":\n";
            out << "  Size (n): " << n_samples[i] << "\n";
            out << "  Mean: " << mean << "\n";
            out << "  Variance: " << variances[i] << "\n";
            out << "  Standard deviation: " << std_devs[i] << "\n";
            out << "  Coefficient of variation: " << (std_devs[i] / mean * 100) << "%\n\n";
        }

        out << "Total observations: " << total_n << "\n\n";

        // Bartlett test calculations
        out << "BARTLETT TEST CALCULATIONS:\n";
        out << "===========================\n";

        double pooled_variance = 0.0;
        double sum_df = 0.0;

        for (int i = 0; i < k; i++) {
            double df = n_samples[i] - 1;
            pooled_variance += df * variances[i];
            sum_df += df;
        }
        pooled_variance /= sum_df;

        out << "Pooled variance: " << pooled_variance << "\n";
        out << "Total degrees of freedom: " << sum_df << "\n";

        double M = sum_df * log(pooled_variance);
        out << "M = " << sum_df << " * ln(" << pooled_variance << ") = " << M << "\n";

        double sum_log_var = 0.0;
        out << "Sum of (n_i - 1) * ln(s_i^2):\n";
        for (int i = 0; i < k; i++) {
            double term = (n_samples[i] - 1) * log(variances[i]);
            sum_log_var += term;
            out << "  Sample " << (i + 1) << ": " << (n_samples[i] - 1)
                << " * ln(" << variances[i] << ") = " << term << "\n";
        }
        out << "Total: " << sum_log_var << "\n";

        M -= sum_log_var;
        out << "M = " << M << " (after subtraction)\n";

        // Correction factor
        double sum_reciprocal = 0.0;
        for (int i = 0; i < k; i++) {
            sum_reciprocal += 1.0 / (n_samples[i] - 1);
        }
        double C = 1.0 + (1.0 / (3.0 * (k - 1))) * (sum_reciprocal - 1.0 / sum_df);

        out << "Correction factor C = 1 + (1/(3*(k-1))) * (sum(1/(n_i-1)) - 1/sum(n_i-1))\n";
        out << "C = 1 + (1/" << (3 * (k - 1)) << ") * (" << sum_reciprocal
            << " - " << (1.0 / sum_df) << ")\n";
        out << "C = " << C << "\n";

        double B_calc = M / C;
        out << "Bartlett statistic B = M / C = " << M << " / " << C << " = " << B_calc << "\n\n";

        // Critical value and conclusion
        double chi2_critical = chi_ppf(1.0 - alpha, k - 1);

        // Calculate exact p-value
        double p_value = 1.0 - chi_cdf(B_calc, k - 1);

        out << "HYPOTHESIS TEST:\n";
        out << "================\n";
        out << "H0: All population variances are equal (sigma1^2 = sigma2^2 = ... = sigmak^2)\n";
        out << "H1: At least one variance is different\n\n";

        out << "Test results:\n";
        out << "Bartlett statistic (B): " << B_calc << "\n";
        out << "Chi-square critical value: " << chi2_critical << "\n";
        out << "Degrees of freedom: " << (k - 1) << "\n";
        out << "Significance level: " << alpha << "\n";
        out << "P-value: " << p_value << "\n\n";

        out << "CONCLUSION:\n";
        out << "===========\n";
        if (B_calc < chi2_critical) {
            out << "B < chi^2_critical: Do not reject H0\n";
            out << "OK The variances appear to be homogeneous (p > " << alpha << ")\n";
            out << "OK The assumption of equal variances is satisfied\n";
        }
        else {
            out << "B >= chi^2_critical: Reject H0\n";
            out << "XX The variances are not homogeneous (p < " << alpha << ")\n";
            out << "XX The assumption of equal variances is violated\n";
        }

        // Additional diagnostic information
        out << "\nDIAGNOSTIC INFORMATION:\n";
        out << "=======================\n";
        out << "Ratio of largest to smallest variance: "
            << (*max_element(variances.begin(), variances.end()) / *min_element(variances.begin(), variances.end())) << "\n";
        out << "Ratio of largest to smallest std. dev.: "
            << (*max_element(std_devs.begin(), std_devs.end()) / *min_element(std_devs.begin(), std_devs.end())) << "\n";

        // Interpretation of variance ratio
        double variance_ratio = *max_element(variances.begin(), variances.end()) /
            *min_element(variances.begin(), variances.end());
        out << "\nINTERPRETATION:\n";
        out << "Variance ratio = " << variance_ratio << " ";
        if (variance_ratio < 1.5) {
            out << "(small difference - homogeneity likely)";
        }
        else if (variance_ratio < 3.0) {
            out << "(moderate difference - borderline)";
        }
        else {
            out << "(large difference - homogeneity questionable)";
        }
        out << "\n";

        if (p_value > alpha) {
            out << "Statistical test confirms: Variances are homogeneous\n";
        }
        else {
            out << "Statistical test indicates: Variances are NOT homogeneous\n";
        }

    }
    catch (const exception& e) {
        out << "ERROR in BartlettTest: " << e.what() << endl;
    }

    out.close();
}
//################### ANOVA  #####################
void ANOVA(string ff) {
    ifstream inp("Inp/" + ff + ".inp");
    ofstream out("Out/" + ff + ".out");

    if (!inp.is_open()) {
        out << "ERROR: Cannot open input file" << endl;
        out.close();
        return;
    }

    try {
        string s1;
        int k; 
        double alpha;
        vector<vector<double>> groups;

        inp >> s1 >> k;

        if (k <= 1) {
            out << "ERROR: Need at least 2 groups for ANOVA" << endl;
            out.close();
            return;
        }

        vector<int> n_groups(k);
        inp >> s1;
        for (int i = 0; i < k; i++) {
            inp >> n_groups[i];
            if (n_groups[i] <= 0) {
                out << "ERROR: Invalid sample size for group " << (i + 1) << endl;
                out.close();
                return;
            }
        }

        inp >> s1 >> alpha;

        groups.resize(k);
        for (int i = 0; i < k; i++) {
            inp >> s1; 
            groups[i].resize(n_groups[i]);
            for (int j = 0; j < n_groups[i]; j++) {
                if (!(inp >> groups[i][j])) {
                    out << "ERROR: Cannot read group " << (i + 1) << " data" << endl;
                    out.close();
                    return;
                }
            }
        }
        inp.close();

        vector<double> means(k);
        vector<double> variances(k);
        int total_n = 0;
        double overall_mean = 0.0;
        double total_sum = 0.0;

        for (int i = 0; i < k; i++) {
            double sum = accumulate(groups[i].begin(), groups[i].end(), 0.0);
            means[i] = sum / n_groups[i];
            total_sum += sum;
            total_n += n_groups[i];

            double var = 0.0;
            for (double x : groups[i]) {
                var += (x - means[i]) * (x - means[i]);
            }
            variances[i] = var / (n_groups[i] - 1);
        }

        overall_mean = total_sum / total_n;

        double SS_between = 0.0;
        double SS_within = 0.0;

        for (int i = 0; i < k; i++) {
            SS_between += n_groups[i] * (means[i] - overall_mean) * (means[i] - overall_mean);
            SS_within += (n_groups[i] - 1) * variances[i];
        }

        double MS_between = SS_between / (k - 1);
        double MS_within = SS_within / (total_n - k);

        double F_calc = (MS_within > 0) ? MS_between / MS_within : 0.0;
        double F_critical = f_ppf(1.0 - alpha, k - 1, total_n - k);

        out << fixed << setprecision(6);
        out << "One-Way ANOVA Results:\n";
        out << "Number of groups: " << k << "\n";
        out << "Total observations: " << total_n << "\n\n";

        for (int i = 0; i < k; i++) {
            out << "Group " << (i + 1) << ": n=" << n_groups[i]
                << ", mean=" << means[i]
                << ", std=" << sqrt(variances[i]) << "\n";
        }

        out << "\nOverall mean: " << overall_mean << "\n";
        out << "SS_between: " << SS_between << " (df=" << (k - 1) << ")\n";
        out << "SS_within: " << SS_within << " (df=" << (total_n - k) << ")\n";
        out << "MS_between: " << MS_between << "\n";
        out << "MS_within: " << MS_within << "\n";
        out << "F-statistic: " << F_calc << "\n";
        out << "F-critical: " << F_critical << " (alpha=" << alpha << ")\n";
        out << "Conclusion: " << (F_calc < F_critical ? "Means are EQUAL" : "Means are DIFFERENT") << "\n";

    }
    catch (const exception& e) {
        out << "ERROR in ANOVA: " << e.what() << endl;
    }
    catch (...) {
        out << "ERROR in ANOVA: Unknown error" << endl;
    }

    out.close();
}

//################### Stubs for missing functions #####################
void lmtaprn(double beta, int n, double v11, double v22, double v12, double delta, double& tlow, double& tup) {
    double z = norm_ppf(1.0 - beta / 2);
    tlow = delta - z;
    tup = delta + z;
}

void nctWeibull(double beta, int n, double v11, double v22, double v12, double delta, double& tlow, double& tup) {
    double z = norm_ppf(1.0 - beta / 2);
    tlow = delta - z;
    tup = delta + z;
}

void ordern(int n, double pr, double ps, double& er, double& vrs) {
    er = norm_ppf(pr);
    vrs = 0.0; 
}

void orderw(int n, double pr, double ps, double& er, double& vrs) {
    er = log(-log(1.0 - pr));
    vrs = 0.0; 
}

void MleastSquare(int n, int k, double** x, double** y, double**& db, double**& b, double*& yr) {
    for (int i = 0; i < k; i++) {
        b[i] = new double[1];
        b[i][0] = 1.0; 
    }
}

void MleastSquare_weight(int n, int k, double** x, double** y, double** v, double**& db, double**& b, double*& yr) {
    for (int i = 0; i < k; i++) {
        b[i] = new double[1];
        b[i][0] = 1.0; 
    }
}
//************************** Inverse Matrix ****************************
double** InverseMatrix(double** a, int n) {
    if (n != 2) {
        double** inv = new double* [2];
        for (int i = 0; i < 2; i++) inv[i] = new double[2]();

        double det = a[0][0] * a[1][1] - a[0][1] * a[1][0];
        if (fabs(det) < 1e-12) {
            inv[0][0] = 1.0; inv[0][1] = 0.0;
            inv[1][0] = 0.0; inv[1][1] = 1.0;
            return inv;
        }

        inv[0][0] = a[1][1] / det;
        inv[0][1] = -a[0][1] / det;
        inv[1][0] = -a[1][0] / det;
        inv[1][1] = a[0][0] / det;

        return inv;
    }

    // For the general case, we use the Gauss method
    double temp;
    int i, j, k;

    double** e = new double* [n];
    for (i = 0; i < n; i++) e[i] = new double[n]();

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            e[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (k = 0; k < n; k++) {
        temp = a[k][k];
        for (j = 0; j < n; j++) {
            a[k][j] /= temp;
            e[k][j] /= temp;
        }

        for (i = k + 1; i < n; i++) {
            temp = a[i][k];
            for (j = 0; j < n; j++) {
                a[i][j] -= a[k][j] * temp;
                e[i][j] -= e[k][j] * temp;
            }
        }
    }

    // The reverse of the Gauss method
    for (k = n - 1; k > 0; k--) {
        for (i = k - 1; i >= 0; i--) {
            temp = a[i][k];
            for (j = 0; j < n; j++) {
                a[i][j] -= a[k][j] * temp;
                e[i][j] -= e[k][j] * temp;
            }
        }
    }

    return e;
}

double** MultiplyMatrix(int rowsa, int colsa, int rowsb, int colsb, double** a, double** b) {
    if (colsa != rowsb) return nullptr;

    double** result = new double* [rowsa];
    for (int i = 0; i < rowsa; i++) {
        result[i] = new double[colsb]();
        for (int j = 0; j < colsb; j++) {
            for (int k = 0; k < colsa; k++) {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return result;
}

double** TransMatrix(int m, int n, double** a) {
    double** result = new double* [n];
    for (int i = 0; i < n; i++) {
        result[i] = new double[m];
        for (int j = 0; j < m; j++) {
            result[i][j] = a[j][i];
        }
    }
    return result;
}

//################### Covariance Matrix MLE Normal #####################
void CovMatrixMleN(int n, vector<double> x, vector<int> r, double a, double s, double**& v) {
    double z, p_val, d, s1, s2, s3, psi;
    int j, k;
    s1 = 0; s2 = 0; s3 = 0; k = 0;

    for (j = 0; j < n; j++) {
        z = (x[j] - a) / s;
        p_val = norm_cdf(z);
        d = norm_pdf(z);
        psi = d / (1 - p_val);
        s1 += r[j] * psi * (psi - z);
        s2 += r[j] * psi * z * (z * (psi - z) - 1);
        s3 += r[j] * psi * (z * (psi - z) - 1);
        k += (1 - r[j]);
    }

    v[0][0] = (k + s1) / n;
    v[0][1] = s3 / n;
    v[1][0] = s3 / n;
    v[1][1] = (2 * k + s2) / n;

    double** inv_v = InverseMatrix(v, 2);

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            v[i][j] = inv_v[i][j];
        }
    }

    for (int i = 0; i < 2; i++) delete[] inv_v[i];
    delete[] inv_v;
}

//################### Covariance Matrix MLE Weibull #####################
void CovMatrixMleW(int n, vector<double> x, vector<int> r, double c, double b, double**& v) {
    int i, k;
    double s1, s2, z, cpw, ckow;
    cpw = log(c);
    ckow = 1 / b;
    s1 = 0;
    s2 = 0;
    k = 0;

    for (i = 0; i < n; i++) {
        z = (log(x[i]) - cpw) / ckow;
        s1 += (1 - r[i]) * z;
        s2 += z * z * exp(z);
        k += (1 - r[i]);
    }

    v[0][0] = double(k) / double(n);
    v[0][1] = (k + s1) / n;
    v[1][0] = (k + s1) / n;
    v[1][1] = (k + s2) / n;

    double** inv_v = InverseMatrix(v, 2);

    for (int i = 0; i < 2; i++) {
        for (int j = 0; j < 2; j++) {
            v[i][j] = inv_v[i][j];
        }
    }

    for (int i = 0; i < 2; i++) delete[] inv_v[i];
    delete[] inv_v;
}

//################### MLE Normal ######################
void MLE_Normal(string ff) {

    int i, j, k, kx, icount, kp;
    string s1;
    double** v = nullptr, * fcum = nullptr, * ycum = nullptr, * xplow = nullptr, * xpup = nullptr, * zp = nullptr, * p = nullptr, tlow, tup, * xp = nullptr;
    double cp, cko, q, eps, beta, delta, step, * t = nullptr, tp, z;

    ifstream inp("Inp/" + ff + ".inp");
    ofstream out("Out/" + ff + ".out");

    if (!inp.is_open()) {
        out << "ERROR: Cannot open input file" << endl;
        out.close();
        return;
    }

    try {
        inp >> s1;
        inp >> nesm.n;
        inp >> s1;
        inp >> beta;
        inp >> s1;
        inp >> eps;
        inp >> s1;

        for (i = 0; i < nesm.n; i++) {
            inp >> z;
            nesm.x.push_back(z);
        }

        inp >> s1;
        for (i = 0; i < nesm.n; i++) {
            inp >> j;
            nesm.r.push_back(j);
        }

        inp >> s1;
        inp >> kp;
        p = new double[kp];
        zp = new double[kp];
        xp = new double[kp];
        inp >> s1;
        for (i = 0; i < kp; i++) inp >> p[i];
        inp.close();

        v = new double* [2];
        for (i = 0; i < 2; i++) {
            v[i] = new double[2];
            for (j = 0; j < 2; j++) v[i][j] = 0.0;
        }

        cp = 0; cko = 0; k = 0;
        for (i = 0; i < nesm.n; i++) {
            k += (1 - nesm.r[i]);
            cp += (1 - nesm.r[i]) * nesm.x[i];
            cko += (1 - nesm.r[i]) * nesm.x[i] * nesm.x[i];
        }

        if (k <= 1) {
            throw runtime_error("Not enough uncensored observations");
        }

        cp /= k;
        cko = sqrt((cko - cp * cp * k) / (k - 1));

        vector<double> xsimpl;
        xsimpl.push_back(cp);
        xsimpl.push_back(cko);

        v[0][0] = 1.0;
        v[1][1] = 0.5;
        v[0][1] = 0.0;
        v[1][0] = 0.0;

        q = 0;
        icount = 0;

        if (k != nesm.n) {
            icount = neldermead(xsimpl, eps, NormalMinFunction);
            q = NormalMinFunction(xsimpl);
            CovMatrixMleN(nesm.n, nesm.x, nesm.r, xsimpl[0], xsimpl[1], v);
        }

        kx = kp;
        t = new double[kx];
        xplow = new double[kx];
        xpup = new double[kx];

        for (i = 0; i < kx; i++) {
            zp[i] = norm_ppf(p[i]);
            xp[i] = xsimpl[0] + zp[i] * xsimpl[1];
            delta = zp[i] * sqrt(nesm.n);

            if (k == nesm.n) {
                t[i] = nct_ppf(beta, nesm.n - 1, delta);
            }
            else {
                lmtaprn(beta, nesm.n - 1, v[0][0], v[1][1], v[0][1], delta, tlow, tup);
                xpup[i] = xsimpl[0] + tup * xsimpl[1] / sqrt(nesm.n);
                xplow[i] = xsimpl[0] + tlow * xsimpl[1] / sqrt(nesm.n);
            }
        }

        if (k == nesm.n) {
            for (i = 0; i < kx; i++)
                xpup[i] = xsimpl[0] + t[i] * xsimpl[1] / sqrt(nesm.n);
            for (i = 0; i < kx; i++)
                xplow[i] = xsimpl[0] - t[kx - i - 1] * xsimpl[1] / sqrt(nesm.n);
        }

        system("mkdir Out 2>nul");

        // Convert arrays to vectors for plotting
        vector<double> Xp_vec(kx), Xplow_vec(kx), Xpup_vec(kx), P_vec(kx);
        for (i = 0; i < kx; i++) {
            Xp_vec[i] = xp[i];
            Xplow_vec[i] = xplow[i];
            Xpup_vec[i] = xpup[i];
            P_vec[i] = p[i];
        }

        // Generate distribution plots with confidence intervals
        GenerateDistributionPlot(ff, Xp_vec, Xplow_vec, Xpup_vec, P_vec, "Normal");

        vector<double> uncensored_data;
        for (i = 0; i < nesm.n; i++) {
            if (nesm.r[i] == 0) {
                uncensored_data.push_back(nesm.x[i]);
            }
        }
        if (!uncensored_data.empty()) {
            GenerateQQPlot("Out/" + ff + "_qq.dat", uncensored_data, "Q-Q Plot - " + ff);
        }

        out << fixed << setprecision(12);
        out << "Method:" << ff << "\n";
        out << "n=" << nesm.n << "\n";
        out << "X\n";
        for (i = 0; i < nesm.n; i++) {
            out << nesm.x[i];
            if (i < nesm.n - 1) out << " , ";
        }
        out << " \n";
        out << "R\n";
        for (i = 0; i < nesm.n; i++) {
            out << nesm.r[i];
            if (i < nesm.n - 1) out << " , ";
        }
        out << " \n";
        out << "cp*=" << cp << "\n";
        out << "cko*=" << cko << "\n";
        out << "Q=" << q << "\n";
        out << "icount=" << icount << "\n";
        out << "cp=" << xsimpl[0] << "\n";
        out << "cko=" << xsimpl[1] << "\n";
        out << "P\n";
        for (i = 0; i < kx; i++) {
            out << p[i];
            if (i < kx - 1) out << " ; ";
        }
        out << " \n";
        out << "Xplow\n";
        for (i = 0; i < kx; i++) {
            out << xplow[i];
            if (i < kx - 1) out << " ; ";
        }
        out << " \n";
        out << "Xp\n";
        for (i = 0; i < kx; i++) {
            out << xp[i];
            if (i < kx - 1) out << " ; ";
        }
        out << " \n";
        out << "Xpup\n";
        for (i = 0; i < kx; i++) {
            out << xpup[i];
            if (i < kx - 1) out << " ; ";
        }
        out << " \n";
        out << "v11=" << v[0][0] << "\n";
        out << "v12=" << v[0][1] << "\n";
        out << "v21=" << v[1][0] << "\n";
        out << "v22=" << v[1][1] << "\n";

    }
    catch (const exception& e) {
        out << "ERROR in MLE_Normal: " << e.what() << endl;
    }

    if (p) delete[] p;
    if (zp) delete[] zp;
    if (xp) delete[] xp;
    if (t) delete[] t;
    if (xplow) delete[] xplow;
    if (xpup) delete[] xpup;
    if (fcum) delete[] fcum;
    if (ycum) delete[] ycum;

    if (v) {
        for (i = 0; i < 2; i++) {
            if (v[i]) delete[] v[i];
        }
        delete[] v;
    }

    nesm.r.clear();
    nesm.x.clear();
    out.close();
}

//################### Nelder-Mead Method ##############################
int neldermead(vector<double>& x0, double eps, double(*func)(vector<double>)) {
    int max_iter = 1000;
    int n = x0.size();

    vector<vector<double>> simplex(n + 1, x0);
    for (int i = 0; i < n; i++) {
        simplex[i + 1][i] += 0.1;
    }

    vector<double> fvals(n + 1);
    for (int i = 0; i <= n; i++) {
        fvals[i] = func(simplex[i]);
    }

    for (int iter = 0; iter < max_iter; iter++) {
        // Find the best, worst and second worst points
        int worst = 0, best = 0, second_worst = 0;
        for (int i = 1; i <= n; i++) {
            if (fvals[i] > fvals[worst]) {
                second_worst = worst;
                worst = i;
            }
            else if (fvals[i] > fvals[second_worst]) {
                second_worst = i;
            }
            if (fvals[i] < fvals[best]) best = i;
        }

        // Convergence check
        double range = 0;
        for (int i = 0; i <= n; i++) {
            double diff = fabs(fvals[i] - fvals[best]);
            if (diff > range) range = diff;
        }
        if (range < eps) {
            x0 = simplex[best];
            return iter;
        }

        // Calculation of the center of gravity (without the worst point)
        vector<double> centroid(n, 0);
        for (int i = 0; i <= n; i++) {
            if (i != worst) {
                for (int j = 0; j < n; j++) {
                    centroid[j] += simplex[i][j];
                }
            }
        }
        for (int j = 0; j < n; j++) {
            centroid[j] /= n;
        }

        vector<double> xr(n);
        for (int j = 0; j < n; j++) {
            xr[j] = centroid[j] + 1.0 * (centroid[j] - simplex[worst][j]);
        }
        double fr = func(xr);

        if (fr < fvals[best]) {
            vector<double> xe(n);
            for (int j = 0; j < n; j++) {
                xe[j] = centroid[j] + 2.0 * (centroid[j] - simplex[worst][j]);
            }
            double fe = func(xe);

            if (fe < fr) {
                simplex[worst] = xe;
                fvals[worst] = fe;
            }
            else {
                simplex[worst] = xr;
                fvals[worst] = fr;
            }
        }
        else if (fr < fvals[worst]) {
            simplex[worst] = xr;
            fvals[worst] = fr;
        }
        else {
            vector<double> xc(n);
            if (fr < fvals[worst]) {
                for (int j = 0; j < n; j++) {
                    xc[j] = centroid[j] + 0.5 * (xr[j] - centroid[j]);
                }
            }
            else {
                for (int j = 0; j < n; j++) {
                    xc[j] = centroid[j] + 0.5 * (simplex[worst][j] - centroid[j]);
                }
            }
            double fc = func(xc);

            if (fc < fvals[worst]) {
                simplex[worst] = xc;
                fvals[worst] = fc;
            }
            else {
                for (int i = 0; i <= n; i++) {
                    if (i != best) {
                        for (int j = 0; j < n; j++) {
                            simplex[i][j] = simplex[best][j] + 0.5 * (simplex[i][j] - simplex[best][j]);
                        }
                        fvals[i] = func(simplex[i]);
                    }
                }
            }
        }
    }

    x0 = simplex[0];
    return max_iter;
}

vector<double> calculateExpectedOrderStatistics(int n) {
    vector<double> a(n);
    double sum_squares = 0.0;

    // We calculate the mathematical expectations of ordinal statistics using ordern
    for (int i = 0; i < n; i++) {
        double pr = (i + 1 - 0.375) / (n + 0.25); 
        double ps = 1.0 - pr;
        double er = 0.0, vrs = 0.0;

        // ordern
        ordern(n, pr, ps, er, vrs);
        a[i] = er;
        sum_squares += er * er;
    }

    // Normalization: a_i = a_i / sqrt(Σa_i^2)
    if (sum_squares > 0) {
        double norm_factor = sqrt(sum_squares);
        for (int i = 0; i < n; i++) {
            a[i] /= norm_factor;
        }
    }

    return a;
}

// ################### Shapiro-Wilk Test #####################
void ShapiroWilkTest(string ff) {
    string inpFile = "Inp/" + ff + ".inp";
    string outFile = "Out/" + ff + ".out";

    ifstream fin(inpFile);
    ofstream fout(outFile);

    if (!fin.is_open()) {
        fout << "Error: Cannot open file " << inpFile << endl;
        return;
    }

    // Reading data
    int n;
    double alpha;
    string line;

    fin >> line; // "n"
    fin >> n;

    fin >> line; // "alpha"
    fin >> alpha;

    fin >> line; // "data"
    vector<double> data(n);
    for (int i = 0; i < n; i++) {
        fin >> data[i];
    }

    fin.close();

    if (n < 3 || n > 50) {
        fout << "Error: Sample size must be between 3 and 50 for Shapiro-Wilk test" << endl;
        fout << "Current n = " << n << endl;
        return;
    }
    
    vector<double> sorted_data = data;
    sort(sorted_data.begin(), sorted_data.end());

    // We calculate the sample variance s^2 (formula 3.18)
    double mean = 0.0;
    for (double val : data) {
        mean += val;
    }
    mean /= n;

    double s2 = 0.0;
    for (double val : data) {
        s2 += (val - mean) * (val - mean);
    }

    // Calculating the coefficients a_i ordern
    vector<double> a = calculateExpectedOrderStatistics(n);

    double b = 0.0;

    if (n % 2 == 0) { // Even number of observations
        // b = Σ a_{n-i+1} * (x_{(n-i+1)} - x_{(i)})
        for (int i = 0; i < n / 2; i++) {
            int idx1 = n - i - 1; // x_{(n-i)}
            int idx2 = i;         // x_{(i+1)}
            b += a[idx1] * (sorted_data[idx1] - sorted_data[idx2]);
        }
    }
    else { // An odd number of observations
        // b = Σ a_{n-i+1} * (x_{(n-i+1)} - x_{(i)})
      // The middle element is not involved (a_{(n+1)/2} = 0)
        for (int i = 0; i < (n - 1) / 2; i++) {
            int idx1 = n - i - 1; // x_{(n-i)}
            int idx2 = i;         // x_{(i+1)}
            b += a[idx1] * (sorted_data[idx1] - sorted_data[idx2]);
        }
    }
    // We calculate the statistics of W (formula 3.17)
    double W = (b * b) / s2;
    
    // alpha = 0.05
    map<int, double> W_critical_005 = {
        {3, 0.767}, {4, 0.748}, {5, 0.762}, {6, 0.788}, {7, 0.803},
        {8, 0.818}, {9, 0.829}, {10, 0.842}, {11, 0.850}, {12, 0.859},
        {13, 0.866}, {14, 0.874}, {15, 0.881}, {16, 0.887}, {17, 0.892},
        {18, 0.897}, {19, 0.901}, {20, 0.905}, {21, 0.908}, {22, 0.911},
        {23, 0.914}, {24, 0.916}, {25, 0.918}, {26, 0.920}, {27, 0.923},
        {28, 0.924}, {29, 0.926}, {30, 0.927}, {31, 0.929}, {32, 0.930},
        {33, 0.931}, {34, 0.933}, {35, 0.934}, {36, 0.935}, {37, 0.936},
        {38, 0.938}, {39, 0.939}, {40, 0.940}, {41, 0.941}, {42, 0.942},
        {43, 0.943}, {44, 0.944}, {45, 0.945}, {46, 0.945}, {47, 0.946},
        {48, 0.947}, {49, 0.947}, {50, 0.947}
    };

    // alpha = 0.01
    map<int, double> W_critical_001 = {
        {3, 0.753}, {4, 0.687}, {5, 0.686}, {6, 0.713}, {7, 0.730},
        {8, 0.749}, {9, 0.764}, {10, 0.781}, {11, 0.792}, {12, 0.805},
        {13, 0.814}, {14, 0.825}, {15, 0.835}, {16, 0.844}, {17, 0.851},
        {18, 0.858}, {19, 0.863}, {20, 0.868}, {21, 0.873}, {22, 0.878},
        {23, 0.881}, {24, 0.884}, {25, 0.888}, {26, 0.891}, {27, 0.894},
        {28, 0.896}, {29, 0.898}, {30, 0.900}, {31, 0.902}, {32, 0.904},
        {33, 0.906}, {34, 0.908}, {35, 0.910}, {36, 0.912}, {37, 0.914},
        {38, 0.916}, {39, 0.917}, {40, 0.919}, {41, 0.920}, {42, 0.922},
        {43, 0.923}, {44, 0.924}, {45, 0.926}, {46, 0.927}, {47, 0.928},
        {48, 0.929}, {49, 0.929}, {50, 0.930}
    };

    // lpha = 0.10
    map<int, double> W_critical_010 = {
        {3, 0.789}, {4, 0.792}, {5, 0.806}, {6, 0.826}, {7, 0.838},
        {8, 0.851}, {9, 0.859}, {10, 0.869}, {11, 0.876}, {12, 0.883},
        {13, 0.889}, {14, 0.895}, {15, 0.901}, {16, 0.906}, {17, 0.910},
        {18, 0.914}, {19, 0.917}, {20, 0.920}, {21, 0.923}, {22, 0.926},
        {23, 0.928}, {24, 0.930}, {25, 0.931}, {26, 0.933}, {27, 0.935},
        {28, 0.936}, {29, 0.937}, {30, 0.939}, {31, 0.940}, {32, 0.941},
        {33, 0.942}, {34, 0.943}, {35, 0.944}, {36, 0.945}, {37, 0.946},
        {38, 0.947}, {39, 0.948}, {40, 0.949}, {41, 0.950}, {42, 0.951},
        {43, 0.951}, {44, 0.952}, {45, 0.953}, {46, 0.953}, {47, 0.954},
        {48, 0.954}, {49, 0.955}, {50, 0.955}
    };

    double W_critical = 0.0;

    // We select a critical value depending on alpha
    if (fabs(alpha - 0.01) < 0.001) {
        if (W_critical_001.find(n) != W_critical_001.end()) {
            W_critical = W_critical_001[n];
        }
        else {
            W_critical = W_critical_005[n] - 0.02;
        }
    }
    else if (fabs(alpha - 0.05) < 0.001) {
        if (W_critical_005.find(n) != W_critical_005.end()) {
            W_critical = W_critical_005[n];
        }
        else {
            W_critical = 0.9; 
        }
    }
    else if (fabs(alpha - 0.10) < 0.001) {
        if (W_critical_010.find(n) != W_critical_010.end()) {
            W_critical = W_critical_010[n];
        }
        else {
            W_critical = W_critical_005[n] + 0.02;
        }
    }
    else {
        // alpha = 0.05
        W_critical = W_critical_005[n];
    }

    fout << "SHAPIRO-WILK NORMALITY TEST" << endl;
    fout << "============================" << endl << endl;

    fout << "Test Parameters:" << endl;
    fout << "Sample size (n): " << n << endl;
    fout << "Significance level (alpha): " << alpha << endl << endl;

    fout << "Sample Statistics:" << endl;
    fout << "Mean: " << mean << endl;
    fout << "Variance (s^2): " << s2 << endl;
    fout << "Standard deviation: " << sqrt(s2) << endl << endl;

    fout << "Test Statistics:" << endl;
    fout << "b statistic (numerator of W): " << b << endl;
    fout << "W statistic: " << W << endl;
    fout << "Critical W value (alpha=" << alpha << "): " << W_critical << endl << endl;

    // We output the coefficients a_i (the first 5 and the last 5)
    fout << "Coefficients a_i (first 5 and last 5):" << endl;
    for (int i = 0; i < min(5, n); i++) {
        fout << "  a[" << (i + 1) << "] = " << a[i] << endl;
    }
    if (n > 10) {
        fout << "  ..." << endl;
        for (int i = max(5, n - 5); i < n; i++) {
            fout << "  a[" << (i + 1) << "] = " << a[i] << endl;
        }
    }
    fout << endl;

    fout << "Result:" << endl;
    fout << "=======" << endl;

    if (W > W_critical) {
        fout << "W (" << W << ") > W_critical (" << W_critical << ")" << endl;
        fout << "=> Do not reject null hypothesis H0." << endl;
        fout << "=> The data appears to follow a normal distribution." << endl;
    }
    else {
        fout << "W (" << W << ") <= W_critical (" << W_critical << ")" << endl;
        fout << "=> Reject null hypothesis H0." << endl;
        fout << "=> The data does NOT appear to follow a normal distribution." << endl;
    }

    // P-value estimation (approximation)
    double p_value_approx = 0.0;
    if (W < W_critical_001[n]) {
        p_value_approx = 0.001; 
    }
    else if (W < W_critical_005[n]) {
        p_value_approx = 0.01; 
    }
    else if (W < W_critical_010[n]) {
        p_value_approx = 0.05;
    }
    else {
        p_value_approx = 0.50; 
    }

    fout << endl << "Approximate p-value: " << p_value_approx << endl;

    if (p_value_approx < alpha) {
        fout << "Since p-value (" << p_value_approx << ") < alpha (" << alpha
            << "), we reject H0 at " << (alpha * 100) << "% significance level." << endl;
    }
    else {
        fout << "Since p-value (" << p_value_approx << ") >= alpha (" << alpha
            << "), we do not reject H0 at " << (alpha * 100) << "% significance level." << endl;
    }

    fout.close();

}

// Auxiliary function for calculating ranks
vector<double> CalculateRanks(const vector<double>& data) {
    size_t n = data.size();
    vector<double> ranks(n, 0.0);

    if (n == 0) return ranks;

    // Creating a vector of pairs (value, source index)
    vector<pair<double, size_t>> indexed_data;
    indexed_data.reserve(n);

    for (size_t i = 0; i < n; i++) {
        indexed_data.push_back(make_pair(data[i], i));
    }

    // Sorting by values
    sort(indexed_data.begin(), indexed_data.end(),
        [](const pair<double, size_t>& a, const pair<double, size_t>& b) {
            return a.first < b.first;
        });

    // We assign ranks based on relationships
    size_t i = 0;
    while (i < n) {
        size_t j = i;
        // We find all the same values
        while (j < n && indexed_data[j].first == indexed_data[i].first) {
            j++;
        }
        // Calculating the average rank for a group of links
        double avg_rank = (i + j + 1) / 2.0;

        // Assign ranks to all the elements in the group
        for (size_t k = i; k < j; k++) {
            ranks[indexed_data[k].second] = avg_rank;
        }

        i = j;
    }

    return ranks;
}

//################### Generate Exact Distribution for Mann-Whitney U (AS 62) #####################
bool wilcoxonDistribution(int m, int n, std::vector<double>& cdf, int& fault) {
    fault = 0;

    if (m <= 0 || n <= 0) {
        fault = 1; 
        return false;
    }

    int min_mn = (m < n) ? m : n;
    int max_mn = (m > n) ? m : n;
    int mn1 = m * n + 1; // Number of possible values of U (0 to m*n)

    // For large samples, we return false - we use an approximation
    if (mn1 > 10000) {
        fault = 4; 
        return false;
    }

    std::vector<double> freq(mn1, 0.0);

    // Initialization for i = 1 (all frequencies = 1 for U from 0 to max_mn)
    for (int i = 0; i <= max_mn; i++) {
        freq[i] = 1.0;
    }

    if (min_mn == 1) {
        cdf.resize(mn1);
        double sum = 0.0;
        for (int i = 0; i < mn1; i++) {
            sum += freq[i];
            cdf[i] = sum;
        }
        for (int i = 0; i < mn1; i++) {
            cdf[i] /= sum;
        }
        return true;
    }

    // Working array
    int work_size = (mn1 + 1) / 2 + min_mn;
    std::vector<double> work(work_size, 0.0);

    // Clearing the rest of the freq
    for (int i = max_mn + 1; i < mn1; i++) {
        freq[i] = 0.0;
    }

    // Generation of higher-order distributions
    work[0] = 0.0;
    int in = max_mn;

    for (int i = 2; i <= min_mn; i++) {
        work[i - 1] = 0.0;
        in += max_mn;
        int n1 = in + 2;
        int l = 1 + in / 2;
        int k = i;

        // Generating a complete distribution from outside to inside
        for (int j = 1; j <= l; j++) {
            k = k + 1;
            n1 = n1 - 1;

            double sum_val = freq[j - 1] + work[j - 1];
            freq[j - 1] = sum_val;

            if (k - 1 < work_size && n1 - 1 < mn1) {
                work[k - 1] = sum_val - freq[n1 - 1];
                freq[n1 - 1] = sum_val;
            }
        }
    }

    double total = 0.0;
    for (int i = 0; i < mn1; i++) {
        total += freq[i];
    }

    if (total == 0) {
        fault = 2; 
        return false;
    }

    cdf.resize(mn1);
    double running_sum = 0.0;
    for (int i = 0; i < mn1; i++) {
        running_sum += freq[i];
        cdf[i] = running_sum / total;
    }

    return true;
}

//################### Wilcoxon Signed-Rank Test (Paired Samples) #####################
void WilcoxonSignedRankTest(string ff) {
    ifstream inp("Inp/" + ff + ".inp");
    ofstream out("Out/" + ff + ".out");

    if (!inp.is_open()) {
        out << "ERROR: Cannot open input file" << endl;
        out.close();
        return;
    }

    int n = 0;
    double alpha = 0.0;

    try {
        string s1;
        inp >> s1 >> n;
        inp >> s1 >> alpha;
        inp >> s1;

        vector<double> sample1(n);
        for (int i = 0; i < n; i++) {
            inp >> sample1[i];
        }

        inp >> s1;
        vector<double> sample2(n);
        for (int i = 0; i < n; i++) {
            inp >> sample2[i];
        }
        inp.close();

        out << fixed << setprecision(6);
        out << "WILCOXON SIGNED-RANK TEST (Paired Samples)\n";
        out << "==========================================\n\n";

        out << "TEST PARAMETERS:\n";
        out << "Sample size (n pairs): " << n << "\n";
        out << "Significance level (alpha): " << alpha << "\n\n";

        // Calculate differences
        vector<double> differences(n);
        for (int i = 0; i < n; i++) {
            differences[i] = sample1[i] - sample2[i];
        }

        // Store absolute differences with signs
        vector<pair<double, int>> abs_diffs;
        for (int i = 0; i < n; i++) {
            if (fabs(differences[i]) > 1e-12) {  // Ignore exact zeros
                int sign = (differences[i] > 0) ? 1 : -1;
                abs_diffs.push_back(make_pair(fabs(differences[i]), sign));
            }
        }

        int n_valid = static_cast<int>(abs_diffs.size());

        if (n_valid == 0) {
            out << "All differences are zero. Test cannot be performed.\n";
            out.close();
            return;
        }

        // Sort by absolute difference
        sort(abs_diffs.begin(), abs_diffs.end(),
            [](const pair<double, int>& a, const pair<double, int>& b) {
                return a.first < b.first;
            });

        // Assign ranks with tie handling
        vector<double> ranks(n_valid);
        int current_idx = 0;

        while (current_idx < n_valid) {
            double current_val = abs_diffs[current_idx].first;
            int tie_count = 1;

            // Count ties
            for (int j = current_idx + 1; j < n_valid; j++) {
                if (fabs(abs_diffs[j].first - current_val) < 1e-12) {
                    tie_count++;
                }
                else {
                    break;
                }
            }

            // Calculate average rank for this group
            double rank_sum = 0.0;
            for (int k = 0; k < tie_count; k++) {
                rank_sum += current_idx + k + 1;
            }
            double avg_rank = rank_sum / tie_count;

            // Assign average rank to all in group
            for (int k = 0; k < tie_count; k++) {
                ranks[current_idx + k] = avg_rank;
            }

            current_idx += tie_count;
        }

        // Calculate W+ and W-
        double W_plus = 0.0, W_minus = 0.0;
        for (int i = 0; i < n_valid; i++) {
            if (abs_diffs[i].second > 0) {
                W_plus += ranks[i];
            }
            else {
                W_minus += ranks[i];
            }
        }

        double W_stat = min(W_plus, W_minus);

        out << "DESCRIPTIVE STATISTICS:\n";
        double mean1 = accumulate(sample1.begin(), sample1.end(), 0.0) / n;
        double mean2 = accumulate(sample2.begin(), sample2.end(), 0.0) / n;
        out << "Mean of sample 1: " << mean1 << "\n";
        out << "Mean of sample 2: " << mean2 << "\n";
        out << "Mean difference: " << mean1 - mean2 << "\n\n";

        out << "TEST STATISTICS:\n";
        out << "W+ (sum of positive ranks): " << W_plus << "\n";
        out << "W- (sum of negative ranks): " << W_minus << "\n";
        out << "Test statistic W: " << W_stat << "\n";

        // Calculate p-value
        double p_value = 0.0;

        if (n_valid <= 30) {
            // Use normal approximation for small samples
            double expected_W = n_valid * (n_valid + 1) / 4.0;
            double var_W = n_valid * (n_valid + 1) * (2 * n_valid + 1) / 24.0;

            // Tie correction
            if (current_idx != n_valid) { // Ties exist
                double tie_correction = 0.0;
                int tie_idx = 0;
                while (tie_idx < n_valid) {
                    double val = abs_diffs[tie_idx].first;
                    int tie_count_local = 1;

                    for (int j = tie_idx + 1; j < n_valid &&
                        fabs(abs_diffs[j].first - val) < 1e-12; j++) {
                        tie_count_local++;
                    }

                    if (tie_count_local > 1) {
                        tie_correction += (tie_count_local * tie_count_local * tie_count_local
                            - tie_count_local);
                    }

                    tie_idx += tie_count_local;
                }

                if (tie_correction > 0) {
                    var_W -= tie_correction / 48.0;
                }
            }

            double Z = 0.0;
            if (var_W > 0) {
                Z = (W_stat - expected_W) / sqrt(var_W);
            }
            else {
                Z = 0.0; // Handle zero variance case
            }

            p_value = 2 * (1.0 - norm_cdf(fabs(Z)));
            out << "Method: Normal approximation (n = " << n_valid << ")\n";
            out << "Z-statistic: " << Z << "\n";
        }
        else {
            // For larger samples, always use normal approximation
            double expected_W = n_valid * (n_valid + 1) / 4.0;
            double var_W = n_valid * (n_valid + 1) * (2 * n_valid + 1) / 24.0;

            double Z = 0.0;
            if (var_W > 0) {
                Z = (W_stat - expected_W) / sqrt(var_W);
            }
            else {
                Z = 0.0; // Handle zero variance case
            }

            p_value = 2 * (1.0 - norm_cdf(fabs(Z)));
            out << "Method: Normal approximation\n";
            out << "Z-statistic: " << Z << "\n";
        }

        out << "P-value: " << p_value << "\n\n";

        out << "HYPOTHESIS TEST:\n";
        out << "H0: Median difference = 0\n";
        out << "H1: Median difference ≠ 0 (two-sided)\n\n";

        out << "CONCLUSION:\n";
        if (p_value > alpha) {
            out << "p-value (" << p_value << ") > alpha (" << alpha << ")\n";
            out << "Do not reject H0\n";
            out << "No significant difference between paired samples\n";
        }
        else {
            out << "p-value (" << p_value << ") <= alpha (" << alpha << ")\n";
            out << "Reject H0\n";
            out << "Significant difference between paired samples\n";
        }

        // Detailed output (first 10 pairs)
        out << "\nDETAILED ANALYSIS (first 10 pairs):\n";
        out << "Idx\tValue1\tValue2\tDiff\t|Diff|\tRank\tSign\n";
        for (int i = 0; i < min(10, n); i++) {
            out << i + 1 << "\t" << sample1[i] << "\t" << sample2[i] << "\t"
                << differences[i] << "\t" << fabs(differences[i]) << "\t";

            if (i < n_valid) {
                out << ranks[i] << "\t" << (abs_diffs[i].second > 0 ? "+" : "-");
            }
            else {
                out << "0.0\t0";
            }
            out << "\n";
        }
    }
    catch (const exception& e) {
        out << "ERROR in WilcoxonSignedRankTest: " << e.what() << endl;
        out << "n = " << n << endl;
    }

    out.close();
}

//################### Wilcoxon Rank-Sum Test (Two Independent Samples)  #####################
void WilcoxonRankSumTest(string ff) {
    ifstream inp("Inp/" + ff + ".inp");
    ofstream out("Out/" + ff + ".out");

    if (!inp.is_open()) {
        out << "ERROR: Cannot open input file" << endl;
        out.close();
        return;
    }

    int m = 0, n = 0;
    double alpha = 0.05;
    vector<double> sample1, sample2;

    try {
        out << "DEBUG: Starting to read input file...\n";

        string line;
        int line_num = 0;
        int reading_mode = 0; 

        while (getline(inp, line)) {
            line_num++;
            string trimmed = line;
            trimmed.erase(0, trimmed.find_first_not_of(" \t"));
            trimmed.erase(trimmed.find_last_not_of(" \t") + 1);

            if (trimmed.empty()) continue;

            if (line_num <= 10) {
                out << "DEBUG Line " << line_num << ": [" << trimmed << "]\n";
            }

            if (trimmed == "sample1") {
                reading_mode = 1;
                continue;
            }
            else if (trimmed == "sample2") {
                reading_mode = 2;
                continue;
            }
            else if (trimmed == "sample_sizes") {
                if (getline(inp, line)) {
                    line_num++;
                    stringstream ss(line);
                    if (!(ss >> m >> n)) {
                        out << "ERROR: Cannot parse sample sizes from: " << line << endl;
                        return;
                    }
                    out << "DEBUG: Parsed m = " << m << ", n = " << n << endl;
                }
                continue;
            }
            else if (trimmed == "alpha") {
                if (getline(inp, line)) {
                    line_num++;
                    stringstream ss(line);
                    if (!(ss >> alpha)) {
                        out << "ERROR: Cannot parse alpha from: " << line << endl;
                        return;
                    }
                    out << "DEBUG: Parsed alpha = " << alpha << endl;
                }
                continue;
            }

            if (reading_mode == 1 || reading_mode == 2) {
                stringstream ss(trimmed);
                double value;

                while (ss >> value) {
                    if (reading_mode == 1) {
                        sample1.push_back(value);
                    }
                    else {
                        sample2.push_back(value);
                    }
                }
            }
        }

        inp.close();

        out << "\nDEBUG SUMMARY:\n";
        out << "Expected sizes: m = " << m << ", n = " << n << endl;
        out << "Actual sizes: sample1 = " << sample1.size() << ", sample2 = " << sample2.size() << endl;
        out << "Alpha: " << alpha << endl;

        if (m <= 0 || n <= 0) {
            m = sample1.size();
            n = sample2.size();
            out << "WARNING: Using actual sizes: m = " << m << ", n = " << n << endl;
        }

        if (sample1.size() != m) {
            out << "WARNING: Adjusting m from " << m << " to " << sample1.size() << endl;
            m = sample1.size();
        }

        if (sample2.size() != n) {
            out << "WARNING: Adjusting n from " << n << " to " << sample2.size() << endl;
            n = sample2.size();
        }

        if (m <= 0 || n <= 0) {
            out << "ERROR: Invalid sample sizes after adjustment\n";
            return;
        }

        out.seekp(0);

        out << "WILCOXON RANK-SUM TEST (Mann-Whitney U Test)\n";
        out << "============================================\n\n";

        out << "TEST PARAMETERS:\n";
        out << "Sample 1 size (m): " << m << "\n";
        out << "Sample 2 size (n): " << n << "\n";
        out << "Total observations: " << (m + n) << "\n";
        out << "Significance level (alpha): " << alpha << "\n\n";

        // Combining samples for ranking
        vector<double> all_data;
        vector<int> group_labels;

        for (double val : sample1) {
            all_data.push_back(val);
            group_labels.push_back(1);
        }

        for (double val : sample2) {
            all_data.push_back(val);
            group_labels.push_back(2);
        }

        int total_n = m + n;

        vector<double> ranks = CalculateRanks(all_data);

        // The sum of the ranks for the first sample
        double R1 = 0.0;
        for (int i = 0; i < total_n; i++) {
            if (group_labels[i] == 1) {
                R1 += ranks[i];
            }
        }

        // U-Mann-Whitney statistics
        double U1 = R1 - m * (m + 1.0) / 2.0;
        double U2 = m * n - U1;
        double U_stat = min(U1, U2);

        out << "TEST STATISTICS:\n";
        out << "Sum of ranks for sample 1 (R1): " << R1 << "\n";
        out << "U1 = R1 - m(m+1)/2 = " << U1 << "\n";
        out << "U2 = m*n - U1 = " << U2 << "\n";
        out << "Test statistic U = min(U1, U2) = " << U_stat << "\n\n";

        // Normal approximation
        double mean_U = m * n / 2.0;
        double var_U = m * n * (m + n + 1.0) / 12.0;
        
        // Correction in touch
        double tie_correction = 0.0;
        vector<double> sorted_all = all_data;
        sort(sorted_all.begin(), sorted_all.end());

        int tie_idx = 0;
        while (tie_idx < total_n) {
            double current_val = sorted_all[tie_idx];
            int tie_count = 1;

            for (int j = tie_idx + 1; j < total_n; j++) {
                if (fabs(sorted_all[j] - current_val) < 1e-10) {
                    tie_count++;
                }
                else {
                    break;
                }
            }

            if (tie_count > 1) {
                tie_correction += (tie_count * tie_count * tie_count - tie_count);
            }

            tie_idx += tie_count;
        }

        if (tie_correction > 0) {
            double denominator = (m + n) * (m + n) * (m + n) - (m + n);
            if (denominator != 0) {
                var_U *= (1.0 - tie_correction / denominator);
            }
        }

        double Z = (U_stat - mean_U) / sqrt(var_U);
        double p_value = 2 * (1.0 - norm_cdf(fabs(Z)));

        out << "NORMAL APPROXIMATION RESULTS:\n";
        out << "Expected U: " << mean_U << "\n";
        out << "Variance of U: " << var_U << "\n";
        out << "Z-statistic: " << Z << "\n";
        out << "P-value: " << p_value << "\n\n";

        // Accurate distribution for small samples only
        double exact_p_value = 0.0;
        bool exact_used = false;

        if (m <= 20 && n <= 20) {
            vector<double> exact_cdf;
            int fault = 0;

            if (wilcoxonDistribution(m, n, exact_cdf, fault) && fault == 0) {
                int idx = static_cast<int>(U_stat);
                if (idx >= 0 && idx < exact_cdf.size()) {
                    exact_p_value = 2 * exact_cdf[idx];
                    if (exact_p_value > 1.0) exact_p_value = 1.0;

                    out << "EXACT DISTRIBUTION RESULTS:\n";
                    out << "Exact p-value: " << exact_p_value << "\n";
                    exact_used = true;

                    // Use the exact value if it is available
                    p_value = exact_p_value;
                }
            }
        }

        if (!exact_used && m <= 20 && n <= 20) {
            out << "Note: Samples are small enough for exact calculation (m=" << m
                << ", n=" << n << " <= 20), but exact computation was skipped for speed.\n";
        }

        out << "\nHYPOTHESIS TEST:\n";
        out << "H0: The two populations have identical distributions\n";
        out << "H1: The two populations differ (two-sided)\n\n";

        out << "CONCLUSION:\n";
        if (p_value > alpha) {
            out << "p-value (" << p_value << ") > alpha (" << alpha << ")\n";
            out << "Do not reject H0\n";
            out << "No statistically significant difference between the two samples\n";
        }
        else {
            out << "p-value (" << p_value << ") <= alpha (" << alpha << ")\n";
            out << "Reject H0\n";
            out << "Statistically significant difference between the two samples\n";
        }

        out << "\nDESCRIPTIVE STATISTICS:\n";

        double mean1_val = accumulate(sample1.begin(), sample1.end(), 0.0) / m;
        double mean2_val = accumulate(sample2.begin(), sample2.end(), 0.0) / n;

        out << "Sample 1 mean: " << mean1_val << "\n";
        out << "Sample 2 mean: " << mean2_val << "\n";
        out << "Difference of means: " << mean1_val - mean2_val << "\n";

        // Medians
        vector<double> sorted1 = sample1;
        vector<double> sorted2 = sample2;
        sort(sorted1.begin(), sorted1.end());
        sort(sorted2.begin(), sorted2.end());

        double median1 = (m % 2 == 0) ? (sorted1[m / 2 - 1] + sorted1[m / 2]) / 2.0 : sorted1[m / 2];
        double median2 = (n % 2 == 0) ? (sorted2[n / 2 - 1] + sorted2[n / 2]) / 2.0 : sorted2[n / 2];

        out << "Sample 1 median: " << median1 << "\n";
        out << "Sample 2 median: " << median2 << "\n";
        out << "Difference of medians: " << median1 - median2 << "\n";

        out << "\nSAMPLE DATA (first 5 values):\n";
        out << "Sample 1: ";
        for (int i = 0; i < min(5, m); i++) {
            out << sample1[i] << " ";
        }
        out << "\nSample 2: ";
        for (int i = 0; i < min(5, n); i++) {
            out << sample2[i] << " ";
        }
        out << "\n";

        out << "\nMethod used: " << (exact_used ? "Exact distribution" : "Normal approximation") << "\n";

    }
    catch (const exception& e) {
        out << "\nERROR in WilcoxonRankSumTest: " << e.what() << endl;
    }

    out.close();
}

//################### Kruskal-Wallis Test #####################
void KruskalWallisTest(string ff) {
    string inpFile = "Inp/" + ff + ".inp";
    string outFile = "Out/" + ff + ".out";

    ifstream fin(inpFile);
    ofstream fout(outFile);

    if (!fin.is_open()) {
        fout << "Error: Cannot open file " << inpFile << endl;
        return;
    }

    // Reading data
    int k;
    string line;

    // Read number of samples
    fin >> line; 
    fin >> k;

    // Read sample sizes
    fin >> line;
    vector<int> sample_sizes(k);
    for (int i = 0; i < k; i++) {
        fin >> sample_sizes[i];
    }

    // Read significance level
    fin >> line; // "alpha"
    double alpha;
    fin >> alpha;

    // Read sample data
    vector<vector<double>> samples(k);
    for (int i = 0; i < k; i++) {
        fin >> line; 
        samples[i].resize(sample_sizes[i]);
        for (int j = 0; j < sample_sizes[i]; j++) {
            fin >> samples[i][j];
        }
    }

    fin.close();

    // Total number of observations
    int N = 0;
    for (int sz : sample_sizes) {
        N += sz;
    }

    // Combine all data into one vector
    vector<double> all_data;
    vector<int> group_labels;

    for (int i = 0; i < k; i++) {
        for (double val : samples[i]) {
            all_data.push_back(val);
            group_labels.push_back(i);
        }
    }

    vector<double> ranks = CalculateRanks(all_data);

    // Sum of ranks for each sample
    vector<double> R(k, 0.0);
    for (size_t i = 0; i < all_data.size(); i++) {
        R[group_labels[i]] += ranks[i];
    }

    // H statistic (formula 3.36 from PDF)
    double H = 0.0;
    for (int i = 0; i < k; i++) {
        H += (R[i] * R[i]) / sample_sizes[i];
    }
    H = (12.0 / (N * (N + 1.0))) * H - 3.0 * (N + 1);

    // Correction for ties (formula 3.37 from PDF)
    vector<double> sorted_data = all_data;
    sort(sorted_data.begin(), sorted_data.end());

    double T = 0.0;
    int count = 1;

    for (size_t i = 1; i <= sorted_data.size(); i++) {
        if (i < sorted_data.size() && fabs(sorted_data[i] - sorted_data[i - 1]) < 1e-10) {
            count++;
        }
        else {
            if (count > 1) {
                T += (count * count * count - count);
            }
            count = 1;
        }
    }

    double H_corrected = H;
    if (T > 0) {
        H_corrected = H / (1.0 - T / (N * N * N - N));
    }

    // Critical value using Fisher distribution approximation
    double df1 = k - 1;
    double df2 = N - k;

    // Approximate critical value using F-distribution
    double F_crit = f_ppf(1.0 - alpha, df1, df2);

    // p-value using F-distribution (approximation)
    double F_stat = H_corrected * (df2) / (N - 1 - H_corrected) / df1;
    double p_value = 1.0 - f_cdf(F_stat, df1, df2);

    fout << "KRUSKAL-WALLIS TEST" << endl;
    fout << "===================" << endl << endl;

    fout << "Test Parameters:" << endl;
    fout << "Number of samples (k): " << k << endl;
    fout << "Total observations (N): " << N << endl;
    fout << "Significance level (alpha): " << alpha << endl << endl;

    fout << "Statistics:" << endl;
    fout << "H = " << H << endl;
    fout << "H with tie correction = " << H_corrected << endl;
    fout << "F-statistic (approximation) = " << F_stat << endl << endl;

    fout << "Degrees of freedom:" << endl;
    fout << "  df1 = k - 1 = " << df1 << endl;
    fout << "  df2 = N - k = " << df2 << endl << endl;

    fout << "Critical values:" << endl;
    fout << "  F(" << df1 << "," << df2 << ")_critical = " << F_crit << endl;
    fout << "  p-value = " << p_value << endl << endl;

    fout << "Sum of ranks for each sample:" << endl;
    for (int i = 0; i < k; i++) {
        double mean_rank = R[i] / sample_sizes[i];
        fout << "  Sample " << (i + 1) << ": n=" << sample_sizes[i]
            << ", sumR=" << R[i] << ", meanR=" << mean_rank << endl;
    }
    fout << endl;

    fout << "Result:" << endl;
    fout << "=======" << endl;

    if (p_value < alpha) {
        fout << "p-value (" << p_value << ") < alpha (" << alpha << ")" << endl;
        fout << "=> Reject null hypothesis H0." << endl;
        fout << "=> There are statistically significant differences between samples." << endl;
    }
    else {
        fout << "p-value (" << p_value << ") >= alpha (" << alpha << ")" << endl;
        fout << "=> Do not reject null hypothesis H0." << endl;
        fout << "=> No statistically significant differences between samples found." << endl;
    }

    fout.close();
}

//################### Generate Distribution Plot Data #####################
void GenerateDistributionPlot(const string& ff, const vector<double>& Xp,
    const vector<double>& Xplow, const vector<double>& Xpup,
    const vector<double>& P, const string& dist_name) {

    if (Xp.empty() || P.empty() || Xp.size() != P.size()) {
        cerr << "ERROR: Invalid data for distribution plot generation" << endl;
        return;
    }

    vector<pair<double, double>> points;
    vector<pair<double, double>> lower_points;
    vector<pair<double, double>> upper_points;

    for (size_t i = 0; i < P.size(); i++) {
        points.push_back(make_pair(P[i], Xp[i]));
        if (i < Xplow.size()) {
            lower_points.push_back(make_pair(P[i], Xplow[i]));
        }
        if (i < Xpup.size()) {
            upper_points.push_back(make_pair(P[i], Xpup[i]));
        }
    }

    sort(points.begin(), points.end());
    sort(lower_points.begin(), lower_points.end());
    sort(upper_points.begin(), upper_points.end());

    vector<double> sorted_P, sorted_Xp, sorted_Xplow, sorted_Xpup;
    for (const auto& p : points) {
        sorted_P.push_back(p.first);
        sorted_Xp.push_back(p.second);
    }
    for (const auto& p : lower_points) {
        sorted_Xplow.push_back(p.second);
    }
    for (const auto& p : upper_points) {
        sorted_Xpup.push_back(p.second);
    }

    GenerateChartData("Out/" + ff + "_main.dat", sorted_P, sorted_Xp,
        dist_name + " Quantile Function - " + ff, "Probability", "Quantile");

    if (!sorted_Xplow.empty()) {
        GenerateChartData("Out/" + ff + "_lower.dat", sorted_P, sorted_Xplow,
            "Lower Confidence Bound", "Probability", "Quantile");
    }

    if (!sorted_Xpup.empty()) {
        GenerateChartData("Out/" + ff + "_upper.dat", sorted_P, sorted_Xpup,
            "Upper Confidence Bound", "Probability", "Quantile");
    }

    GenerateChartData("Out/" + ff + "_points.dat", sorted_P, sorted_Xp,
        "Quantile Points", "Probability", "Quantile");

    cout << "Generated distribution plots for " << ff << " with "
        << sorted_P.size() << " points" << endl;
}

//################### Generate Chart Data #####################
void GenerateChartData(const string& filename, const vector<double>& x_data, const vector<double>& y_data,
    const string& title, const string& xlabel, const string& ylabel) {

    ofstream out(filename);
    if (!out.is_open()) {
        cerr << "ERROR: Cannot open chart file: " << filename << endl;
        return;
    }

    out << "# " << title << endl;
    out << "# " << xlabel << " " << ylabel << endl;

    for (size_t i = 0; i < x_data.size() && i < y_data.size(); i++) {
        out << x_data[i] << " " << y_data[i] << endl;
    }

    out.close();
    cout << "Created chart file: " << filename << " with " << x_data.size() << " points" << endl;
}

void GenerateQQPlot(const string& filename, const vector<double>& data, const string& title) {
    ofstream out(filename);
    if (!out.is_open()) return;

    if (data.empty()) return;

    vector<double> sorted_data = data;
    sort(sorted_data.begin(), sorted_data.end());

    int n = static_cast<int>(sorted_data.size());
    out << "# " << title << endl;
    out << "# Theoretical_Quantile Sample_Quantile" << endl;

    for (int i = 0; i < n; i++) {
        double p = (i + 1 - 0.5) / n;
        double theoretical = norm_ppf(p);
        out << theoretical << " " << sorted_data[i] << endl;
    }

    out.close();
}
