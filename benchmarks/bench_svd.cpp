

#include <opencv2/opencv.hpp>
#include <Eigen/Dense>
#include "squareMatrix.h"

#include <filesystem>
#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>

namespace fs = std::filesystem;
using Clock    = std::chrono::high_resolution_clock;
using Duration = std::chrono::duration<double>;

static double frobenius(const linalg::squareMatrix<double>& A,
                        const linalg::squareMatrix<double>& B) {
    size_t n = A.rows();
    double sum = 0.0;
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            double d = A(i,j) - B(i,j);
            sum += d * d;
        }
    }
    return std::sqrt(sum);
}

int main(int argc, char** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <input_dir> <output_dir> <ranks_csv>\n";
        return 1;
    }
    fs::path in_dir  = argv[1];
    fs::path out_dir = argv[2];
    std::string ranks_str = argv[3];

    std::vector<int> ranks;
    {
        std::stringstream ss(ranks_str);
        for (std::string tok; std::getline(ss, tok, ','); ) {
            ranks.push_back(std::stoi(tok));
        }
    }

    fs::create_directories(out_dir / "my");
    fs::create_directories(out_dir / "eigen");

    std::ofstream csv(out_dir / "results.csv");
    csv << "filename,my_svd_time,eigen_svd_time,k,my_recon_time,eigen_recon_time,my_err,eigen_err\n";

    for (auto const& entry : fs::directory_iterator(in_dir)) {
        if (!entry.is_regular_file()) continue;
        auto path = entry.path();
        auto ext  = path.extension().string();
        if (ext != ".png" && ext != ".jpg" && ext != ".jpeg") continue;

        cv::Mat I = cv::imread(path.string(), cv::IMREAD_GRAYSCALE);
        if (I.empty()) continue;
        int m = I.rows, n = I.cols;
        if (m != n) {
            std::cerr << "Skipping non-square: " << path << "\n";
            continue;
        }

        cv::Mat Id;
        I.convertTo(Id, CV_64F, 1.0 / 255.0);

        size_t a = m;
        linalg::vector<double> data(a * a);
        linalg::squareMatrix<double> A(a, data);

        auto t0 = Clock::now();
        auto svd_res = A.svd();
        auto t1 = Clock::now();
        double my_svd_time = Duration(t1 - t0).count();

        auto U  = std::get<0>(svd_res);
        auto S  = std::get<1>(svd_res);
        auto VT = std::get<2>(svd_res);

        Eigen::MatrixXd E(a, a);
        for (int i = 0; i < a; ++i)
            for (int j = 0; j < a; ++j)
                E(i,j) = Id.at<double>(i,j);

        auto t0e = Clock::now();
        Eigen::JacobiSVD<Eigen::MatrixXd> esvd(
                E, Eigen::ComputeThinU | Eigen::ComputeThinV
        );
        auto t1e = Clock::now();
        double eigen_svd_time = Duration(t1e - t0e).count();

        for (int k : ranks) {
            auto Sk = S;
            for (int i = k; i < (int)a; ++i)
                Sk(i,i) = 0;

            auto t2 = Clock::now();
            auto Arec_ten = U.multiply(Sk).multiply(VT);

            auto t3 = Clock::now();
            double my_recon_time = Duration(t3 - t2).count();
            linalg::squareMatrix<double> Arec_my(a, Arec_ten.get_data());
            double my_err        = frobenius(A, Arec_my);

            cv::Mat Imy(a, a, CV_64F);
            for (int i = 0; i < (int)a; ++i)
                for (int j = 0; j < (int)a; ++j)
                    Imy.at<double>(i,j) = Arec_my(i,j);
            Imy.convertTo(Imy, CV_8U, 255.0);
            Imy = Imy.clone();
            cv::imwrite((out_dir / "my" /
                         (path.stem().string() + "_k" + std::to_string(k) + ext)
                        ).string(), Imy);

            auto t2e = Clock::now();
            Eigen::MatrixXd Ue = esvd.matrixU().leftCols(k);
            Eigen::MatrixXd Se = esvd.singularValues().head(k).asDiagonal();
            Eigen::MatrixXd Ve = esvd.matrixV().leftCols(k).transpose();
            Eigen::MatrixXd Erec = Ue * Se * Ve;
            auto t3e = Clock::now();
            double eigen_recon_time = Duration(t3e - t2e).count();
            double eigen_err        = (E - Erec).norm();

            cv::Mat Iei(a, a, CV_64F);
            for (int i = 0; i < (int)a; ++i)
                for (int j = 0; j < (int)a; ++j)
                    Iei.at<double>(i,j) = Erec(i,j);
            Iei.convertTo(Iei, CV_8U, 255.0);
            cv::imwrite((out_dir / "eigen" /
                         (path.stem().string() + "_k" + std::to_string(k) + ext)
                        ).string(), Iei);

            csv << path.filename().string() << ","
                << my_svd_time     << ","
                << eigen_svd_time  << ","
                << k               << ","
                << my_recon_time   << ","
                << eigen_recon_time<< ","
                << my_err          << ","
                << eigen_err       << "\n";
        }
    }

    csv.close();
    return 0;
}
