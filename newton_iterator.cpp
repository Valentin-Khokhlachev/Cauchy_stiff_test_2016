#include "newton_iterator.h"

using namespace Newton;

Newton_iterator::Newton_iterator()
{
    srand(time(nullptr));

    std::vector<std::string> Newton_config_str = parse_config_file("../../../vpa_stiff_test_2016_arg_arc/config/config.txt", ' ', "Newton_config");

    if (DEBUG_MOD) {
        for (int i = 0; i < Newton_config_str.size(); ++i) {
            qDebug() << Newton_config_str[i].c_str();
        }
        qDebug() << Newton_config_str[Newton_config_str.size()-1].c_str();
    }

    std::istringstream(Newton_config_str[0]) >> precForPrint;
    eps = mpComplex_t{Newton_config_str[1], "0.0"};
    theta = mpComplex_t{Newton_config_str[2], "0.0"};
    std::istringstream(Newton_config_str[3]) >> maxIterations;
    alpha = mpComplex_t{Newton_config_str[4], "0.0"};

}

void Newton_iterator::setPrecForPrint(int prc)
{
    precForPrint = prc;
}

int Newton_iterator::getPrecForPrint()
{
    return precForPrint;
}

void Newton_iterator::setEps(string realPart, string imagPart)
{
    this->eps = mpComplex_t{realPart, imagPart};
}

mpComplex_t Newton_iterator::getEps()
{
    return this->eps;
}

void Newton_iterator::setTheta(string realPart, string imagPart)
{
    this->theta = mpComplex_t{realPart, imagPart};
}

mpComplex_t Newton_iterator::getTheta()
{
    return this->theta;
}

void Newton_iterator::setMaxIterations(int maxIt)
{
    this->maxIterations = maxIt;
}

int Newton_iterator::getMaxIterations()
{
    return this->maxIterations;
}

void Newton_iterator::setAlpha(string realPart, string imagPart)
{
    this->alpha = mpComplex_t{realPart, imagPart};
}

mpComplex_t Newton_iterator::getAlpha()
{
    return this->alpha;
}

VectorX<mpComplex_t> Newton_iterator::Newton_method(function<VectorX<mpComplex_t> (VectorX<mpComplex_t>)> f,
                                                     function<MatrixX<mpComplex_t> (VectorX<mpComplex_t>)> J_f,
                                                     const VectorX<mpComplex_t> &v_prev)
{
    int N = v_prev.size();
    VectorX<mpComplex_t> x_next = v_prev;
    VectorX<mpComplex_t> x_prev = v_prev;
    VectorX<mpComplex_t> x_pprev = v_prev;
    VectorX<mpComplex_t> f_prev(N);
    f_prev.setZero();
    VectorX<mpComplex_t> f_prev_xi(N);
    f_prev_xi.setZero();
    MatrixX<mpComplex_t> J_f_prev(N, N);
    J_f_prev.setZero();
    VectorX<mpComplex_t> xi(N);
    xi.setZero();

    mpComplex_t xiNorm {"0.0", "0.0"};
    mpComplex_t phi0 {"0.0", "0.0"};
    mpComplex_t phi1 {"0.0", "0.0"};
    mpComplex_t tau {"0.0", "0.0"};
    mpComplex_t q {"0.0", "0.0"};
    mpComplex_t left {"0.0", "0.0"};
    mpComplex_t right {"0.0", "0.0"};

    int iterator = 0;
    bool flag1 = true;
    bool flag2 = true;

    while(flag1 && flag2){

        qDebug() << "===+===+===+===+===+===+===+===+===+===+===+===";
        qDebug() << "iteration = " << iterator;

        x_pprev = x_prev;
        x_prev = x_next;
        f_prev = f(x_prev);
        J_f_prev = J_f(x_prev);
        xi = J_f_prev.colPivHouseholderQr().solve(f_prev);
        xiNorm = sqrt(static_cast<mpComplex_t>(xi.adjoint() * xi));
        f_prev_xi = f(x_prev + xi);

        qDebug() << "======";
        qDebug() << "xi = " << xiNorm.str(precForPrint).c_str();

        phi0 = static_cast<mpComplex_t>(f_prev.adjoint() * f_prev);
        phi1 = static_cast<mpComplex_t>(f_prev_xi.adjoint() * f_prev_xi);

        qDebug() << "======";
        qDebug() << "phi0 = " << phi0.str(precForPrint).c_str();
        qDebug() << "phi1 = " << phi1.str(precForPrint).c_str();


        tau = (phi0 + phi1 * theta)/(phi0 + phi1);
        x_next = x_prev - tau * xi;

        qDebug() << "======";
        qDebug() << "tau = " << tau.str(precForPrint).c_str();

        if(iterator > 0){
            q = sqrt(abs(static_cast<mpComplex_t>(
                    (x_next.adjoint() - x_prev.adjoint())*(x_next - x_prev))))/
                sqrt(abs(static_cast<mpComplex_t>(
                    (x_prev.adjoint() - x_pprev.adjoint())*(x_prev - x_pprev))));

            qDebug() << "======";
            qDebug() << "q = " << q.str(precForPrint).c_str();

            left = abs(q / (mpComplex_t{"1.0", "0.0"} - q) *
                       sqrt(abs(static_cast<mpComplex_t>(
                           (x_next.adjoint() - x_prev.adjoint()) * (x_next - x_prev))))
                       );

            qDebug() << "======";
            qDebug() << "left = " << left.str(precForPrint).c_str();

            // неправильно, исправить на одновременное комплексное сопряжение и транспонирование
            right = abs(eps * sqrt(abs(static_cast<mpComplex_t>(x_next.adjoint() * x_next))));

            qDebug() << "======";
            qDebug() << "right = " << right.str(precForPrint).c_str();

            if(left.real() > right.real()){
                flag1 = true;
            }else{
                flag1 = false;
            }
        }else{
            flag1 = true;
        }

        if(iterator < maxIterations){
            flag2 = true;
        }else{
            flag2 = false;
        }

        if(iterator == maxIterations){
            qDebug() << "======";
            qDebug() << "maxIterations reached!";
        }

        qDebug() << "===+===+===+===+===+===+===+===+===+===+===+===";

        ++iterator;
    }
    return x_next;
}
