#ifndef DATAGENERATOR_H
#define DATAGENERATOR_H

#include <string>
#include <vector>
#include <random>
#include <ctime>
#include <algorithm>
#include <numeric>
#include "stat.h"

using namespace std;

class DataGenerator {
private:
    static default_random_engine& getGenerator() {
        static default_random_engine generator(static_cast<unsigned int>(time(0)));
        return generator;
    }

public:
    static double randomDouble(double min, double max) {
        uniform_real_distribution<double> distribution(min, max);
        return distribution(getGenerator());
    }

    static int randomInt(int min, int max) {
        uniform_int_distribution<int> distribution(min, max);
        return distribution(getGenerator());
    }

    // Generating data for MLE Normal
    static string generateMLE_Normal() {
        int n = 25;
        double a = 10.0, sigma = 2.0;
        string content = "sample_size\n" + to_string(n) +
            "\nbeta\n0.05\neps\n0.001\nx\n";

        for (int i = 0; i < n; i++) {
            double z = randomDouble(0.001, 0.999);
            double znorm = norm_ppf(z);
            double value = a + znorm * sigma;
            content += to_string(value) + " ";
        }

        content += "\nr\n";
        for (int i = 0; i < n; i++) {
            content += to_string(randomInt(0, 1)) + " ";
        }

        content += "\nkp\n15\np\n";
        vector<double> p_values = { 0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.975, 0.99, 0.995 };
        for (double p_val : p_values) {
            content += to_string(p_val) + " ";
        }

        return content;
    }

    // Generating data for MLE Weibull
    static string generateMLE_Weibull() {
        int n = 25;
        double b = 2.0, c = 1.5;

        string content = "sample_size\n" + to_string(n) +
            "\nbeta\n0.05\neps\n0.001\nx\n";

        for (int i = 0; i < n; i++) {
            double z = randomDouble(0.001, 0.999);
            double value = c * pow(-log(1.0 - z), 1.0 / b);
            content += to_string(value) + " ";
        }

        content += "\nr\n";
        for (int i = 0; i < n; i++) {
            content += "0 ";
        }

        content += "\nkp\n15\np\n";
        vector<double> p_values = { 0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.975, 0.99, 0.995 };
        for (double p_val : p_values) {
            content += to_string(p_val) + " ";
        }

        return content;
    }

    // Data generation for Shapiro-Wilk Test
    static string generateShapiroWilkTest() {
        int n = 25;
        double a = 10.0, sigma = 2.0;

        string content = "n\n" + to_string(n) +
            "\nalpha\n0.05\ndata\n";

        for (int i = 0; i < n; i++) {
            double z = randomDouble(0.001, 0.999);
            double znorm = norm_ppf(z);
            double value = a + znorm * sigma;
            content += to_string(value) + " ";
        }

        return content;
    }

    // Generating data for Wilcoxon Test
    static string generateWilcoxonSignedRankTest() {
        int sample_size = 25;
        double mean1 = 10.0, std_dev1 = 2.0;
        double mean2 = 10.5, std_dev2 = 2.0; 

        string content = "n\n" + to_string(sample_size) + "\nalpha\n0.05\ndata1\n";

        // First selection (before)
        for (int i = 0; i < sample_size; i++) {
            double uniform_random = randomDouble(0.001, 0.999);
            double z_score = norm_ppf(uniform_random);
            double value = mean1 + z_score * std_dev1;
            content += to_string(value) + " ";
        }
        content += "\ndata2\n";

        // The second selection (after)
        for (int i = 0; i < sample_size; i++) {
            double uniform_random = randomDouble(0.001, 0.999);
            double z_score = norm_ppf(uniform_random);
            double value = mean2 + z_score * std_dev2;
            content += to_string(value) + " ";
        }
        return content;
    }

    // For the two-choice criterion (Mann-Whitney)
    static string generateWilcoxonRankSumTest() {
        int m = 15; 
        int n = 20; 
        double mean1 = 10.0, std_dev1 = 2.0;
        double mean2 = 12.0, std_dev2 = 2.0;

        string content = "sample_sizes\n" + to_string(m) + " " + to_string(n) +
            "\nalpha\n0.05\nsample1\n";

        // First selection 
        for (int i = 0; i < m; i++) {
            double uniform_random = randomDouble(0.001, 0.999);
            double z_score = norm_ppf(uniform_random);
            double value = mean1 + z_score * std_dev1;
            content += to_string(value) + " ";
        }
        content += "\nsample2\n";

        // The second selection
        for (int i = 0; i < n; i++) {
            double uniform_random = randomDouble(0.001, 0.999);
            double z_score = norm_ppf(uniform_random);
            double value = mean2 + z_score * std_dev2;
            content += to_string(value) + " ";
        }

        return content;
    }

    // Data generation for MLS Normal and Weibull
    static string generateMLS_Normal() {
        return generateMLS_Base(0); 
    }

    static string generateMLS_Weibull() {
        return generateMLS_Base(1); 
    }

private:
    static string generateMLS_Base(int dist_type) {
        int n = 25;
        double a = 10.0, sigma = 2.0;
        double b = 2.0, c = 1.5;

        string content = "n\n" + to_string(n) + "\nalpha\n0.05\ndata\n";

        for (int i = 0; i < n; i++) {
            double value;
            if (dist_type == 0) { // Normal
                double z = randomDouble(0.001, 0.999);
                double znorm = norm_ppf(z);
                value = a + znorm * sigma;
            }
            else { // Weibull
                double z = randomDouble(0.001, 0.999);
                value = c * pow(-log(1.0 - z), 1.0 / b);
            }
            content += to_string(value) + " ";
        }

        // For MLS, we need r tags (all uncensored)
        content += "\nr\n";
        for (int i = 0; i < n; i++) {
            content += "0 ";
        }

        // And quantiles
        content += "\nkp\n15\np\n";
        vector<double> p_values = { 0.005, 0.01, 0.025, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.8, 0.9, 0.95, 0.975, 0.99, 0.995 };
        for (double p_val : p_values) {
            content += to_string(p_val) + " ";
        }

        return content;
    }

public:
    // Generating data for Grubbs Test
    static string generateGrubbsTest() {
        int n = 25;
        double a = 10.0, sigma = 2.0;

        string content = "n\n" + to_string(n) + "\nalpha\n0.05\ndata\n";

        for (int i = 0; i < n - 1; i++) {
            double z = randomDouble(0.001, 0.999);
            double znorm = norm_ppf(z);
            double value = a + znorm * sigma;
            content += to_string(value) + " ";
        }

        double outlier = a + 5 * sigma;
        content += to_string(outlier);

        return content;
    }

    // Generating data for the Student Test
    static string generateStudentTest() {
        int n1 = 20, n2 = 20;
        double a1 = 10.0, sigma1 = 2.0;
        double a2 = 12.0, sigma2 = 2.0;

        string content = "n1 n2\n" + to_string(n1) + " " + to_string(n2) +
            "\nalpha\n0.05\nsample1\n";

        for (int i = 0; i < n1; i++) {
            double z = randomDouble(0.001, 0.999);
            double znorm = norm_ppf(z);
            double value = a1 + znorm * sigma1;
            content += to_string(value) + " ";
        }

        content += "\nsample2\n";
        for (int i = 0; i < n2; i++) {
            double z = randomDouble(0.001, 0.999);
            double znorm = norm_ppf(z);
            double value = a2 + znorm * sigma2;
            content += to_string(value) + " ";
        }

        return content;
    }

    // Generating data for the Bartlett Test
    static string generateBartlettTest() {
        int k = 3; 
        int n1 = 20, n2 = 20, n3 = 20;

        string content = "k\n" + to_string(k) +
            "\nn_samples\n" + to_string(n1) + " " + to_string(n2) + " " + to_string(n3) +
            "\nalpha\n0.05\nsample1\n";

        // Generating three samples
        for (int i = 0; i < n1; i++) {
            double value = 10.0 + randomDouble(-5.0, 5.0);
            content += to_string(value) + " ";
        }

        content += "\nsample2\n";
        for (int i = 0; i < n2; i++) {
            double value = 12.0 + randomDouble(-5.0, 5.0);
            content += to_string(value) + " ";
        }

        content += "\nsample3\n";
        for (int i = 0; i < n3; i++) {
            double value = 11.0 + randomDouble(-5.0, 5.0);
            content += to_string(value) + " ";
        }

        return content;
    }

    // Generating data for ANOVA
    static string generateANOVA() {
        int k = 3; 
        int n1 = 15, n2 = 15, n3 = 15;

        string content = "k\n" + to_string(k) +
            "\nn_groups\n" + to_string(n1) + " " + to_string(n2) + " " + to_string(n3) +
            "\nalpha\n0.05\ngroup1\n";

        for (int i = 0; i < n1; i++) {
            double value = 10.0 + randomDouble(-3.0, 3.0);
            content += to_string(value) + " ";
        }

        content += "\ngroup2\n";
        for (int i = 0; i < n2; i++) {
            double value = 12.0 + randomDouble(-3.0, 3.0);
            content += to_string(value) + " ";
        }

        content += "\ngroup3\n";
        for (int i = 0; i < n3; i++) {
            double value = 14.0 + randomDouble(-3.0, 3.0);
            content += to_string(value) + " ";
        }

        return content;
    }

    // Generating data for generateKruskalWallis
    static string generateKruskalWallisTest() {
        int k = 5; 
        vector<int> sample_sizes = { 25, 30, 28, 32, 27 }; 

        string content = "k\n" + to_string(k) +
            "\nsample_sizes\n";

        for (size_t i = 0; i < sample_sizes.size(); i++) {
            content += to_string(sample_sizes[i]);
            if (i < sample_sizes.size() - 1) content += " ";
        }

        content += "\nalpha\n0.05\n";

        // Generating data for each sample with different parameters
        vector<double> means = { 50.0, 55.0, 60.0, 65.0, 70.0 }; // different averages
        vector<double> sigmas = { 5.0, 6.0, 7.0, 8.0, 9.0 }; // different standard deviations

        for (int group = 0; group < k; group++) {
            content += "sample\n"; 

            for (int i = 0; i < sample_sizes[group]; i++) {
                double z = randomDouble(0.001, 0.999);
                double znorm = norm_ppf(z);
                double value = means[group] + znorm * sigmas[group];
                value += randomDouble(-1.0, 1.0);

                content += to_string(value);
                if (i < sample_sizes[group] - 1) content += " ";
            }

            if (group < k - 1) content += "\n";
        }

        return content;
    }
};


#endif