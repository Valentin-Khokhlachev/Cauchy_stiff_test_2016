#include "ivp_cauchy_solver.h"

using namespace mpTypes;
using namespace Cauchy;

ivpCauchySolver::ivpCauchySolver()
{
    // Подготовка к генерации псевдослучайных чисел
    srand(time(nullptr));

    // Ввод входных параметров
    std::vector<std::string> Cauchy_config_str = parse_config_file("../../../vpa_stiff_test_2016_arg_arc/config/config.txt", ' ', "Cauchy_config");

    if (DEBUG_MOD) {
        for (int i = 0; i < Cauchy_config_str.size(); ++i) {
            qDebug() << Cauchy_config_str[i].c_str();
        }
        qDebug() << Cauchy_config_str[Cauchy_config_str.size()-1].c_str();
    }

    std::istringstream(Cauchy_config_str[0]) >> N;
    l_0 = mpComplex_t{Cauchy_config_str[1], "0.0"};
    L = mpComplex_t{Cauchy_config_str[2], "0.0"};
    t_0 = mpComplex_t{Cauchy_config_str[3], "0.0"};;
    alpha = mpComplex_t{Cauchy_config_str[4], "0.0"};;
    pA = mpComplex_t{Cauchy_config_str[5], "0.0"};;
    pLambda0 = mpComplex_t{Cauchy_config_str[6], "0.0"};;
    filePath = Cauchy_config_str[7];

    // Подготовка к выводу результатов
    // Если папки для вывода результатов нет, то создаём её
    if(!std::filesystem::is_directory(filePath)){
        std::filesystem::create_directory(filePath);
        qDebug() << "Target folder created!";
    }else{
        qDebug() << "Target folder already exists!";
    }
    if(std::filesystem::is_directory(filePath)){
        // Создаём экземпляр класса csv-файл
        // csvfile = make_unique<csvstream>(file);
        // Открываем csv файл для вывода
        ofs.open(file);
        if(ofs.is_open()){
            qDebug() << '"' << file.c_str() << '"' << "is opened!";
        }else{
            qDebug() << '"' << file.c_str() << "\":" << "can't open!";
        }
    }else{
        qDebug() << "Folder doesn't exist!";
    }

    // Выделение памяти под Ньютоновский модуль
    iter = make_unique<Newton::Newton_iterator>();

    // Инициализация сеток
    phi.resize(N);
    phi = Eigen::VectorX<mpComplex_t>::LinSpaced(N,
                                                  mpComplex_t("0.0"),
                                                  mpComplex_t("2.0") * PI_C);

    l.resize(N);
    l = Eigen::VectorX<mpComplex_t>::LinSpaced(N, l_0, L);

    z.resize(N);
    z = contour(l_0, L, rho, phi);

    h_z.resize(N);
    h_z = hContour(rho, phi, h_phi);

    // Инициализация приближённого решения
    u.resize(2 * N);
    u.setZero();

    // Инициализация начальных условий
    uInit.resize(2);
    uInit[0] = mpFloat_t("-2.0") * (pLambda0 * sin(t_0) * pA * pA /
                                  (mpFloat_t("1.0") + sqrt(mpFloat_t("1.0") +
                                                           mpFloat_t("4.0") * pLambda0 * pLambda0 *
                                                               sin(t_0) * sin(t_0) * pA * pA)));
    uInit[1] = t_0;

    // Инициализация сеточного аналитического продолжения
    uWave.resize(2 * N);
    uWave.setZero();

    // Инициализация начального приближения для сеточного аналитического продолжения
    uWaveInit.resize(2 * N);
    uWaveInit.setZero();

    // Инициализация точного решения
    uExact.resize(2 * N);
    uExact.setZero();
}

ivpCauchySolver::~ivpCauchySolver()
{
    if(ofs.is_open()){
        ofs.close();
        qDebug() << '"' << file.c_str() << '"' << "is closed!";
    }
}

void ivpCauchySolver::solve()
{
    // MAIN PART
    // генерируем начальное приближение для сеточного аналитического продолжения
    this->NewtonInit();
    // лямбда-функция для вызова F
    auto lF = [this](VectorX<mpComplex_t> x)->VectorX<mpComplex_t>{return this->F(x);};
    // лямбда-функция для вызова JF
    auto lJF = [this](VectorX<mpComplex_t> x)->MatrixX<mpComplex_t>{return this->JF(x);};
    // запускаем итерационный процесс
    uWave = this->iter->Newton_method(lF, lJF, uWaveInit);
    // пересчитываем сеточное аналитическое продолжение в искомую сеточную функцию
    this->reCalc();
    // считаем точное решение
    this->exactSolution(u.segment(N, N));
    // считаем относительную погрешность
    mpComplex_t err {"0.0", "0.0"};
    err = this->relErr();
    qDebug() << "Relative error =" << err.str().c_str();
}

VectorX<mpComplex_t> ivpCauchySolver::contour(const mpComplex_t t_s,
                                                 const mpComplex_t t_f,
                                                 const mpComplex_t rho,
                                                 const Ref<VectorX<mpComplex_t>> phi)
{
    Eigen::VectorX<mpComplex_t> res(N);
    for (int i = 0; i < N; ++i){
        res[i] =(t_f + t_s) / 2.0 + exp(I * phi[i]) * rho;
    }
    return res;
}

VectorX<mpComplex_t> ivpCauchySolver::hContour(const mpComplex_t rho,
                                                   const VectorX<mpComplex_t> phi,
                                                   const mpComplex_t h_phi)
{
    Eigen::VectorX<mpComplex_t> res(N);
    for (int i = 0; i < N; ++i){
        res[i] =  I * h_phi * exp(I * phi[i]) * rho;
    }
    return res;
}

string ivpCauchySolver::mstr(const mpComplex_t z, const int pr)
{
    string res = "\0";
    string temp = "\0";
    temp = z.real().str(pr);
    res += temp;
    temp = z.imag().str(pr);
    if(temp[0] == '-'){
        res += temp;
    }else{
        if(temp[0] == '+'){
            res += temp;
        }else{
            res += "+";
            res += temp;
        }
    }
    res += "i";
    return res;
}

void ivpCauchySolver::NewtonInit()
{
    VectorX<mpComplex_t> res(2 * N);
    res.setZero();
    MatrixX<mpComplex_t> temp_matrix(N, N);
    temp_matrix.setZero();

    for(int j = 0; j < N; ++j){
        for(int i = 0; i < N; ++i){
            temp_matrix(j, i) = h_z[i]/(z[i] - l[j]);
        }
    }

    VectorX<mpComplex_t> temp_vec_1(N);
    temp_vec_1.setZero();

    // первый этап
    for(int j = 0; j < N; ++j){
        temp_vec_1[j] = mpFloat_t("2.0") * PI * I *
                        (uInit[0] + rightHand0(uInit) * (l[j] - l_0));
    }

    VectorX<mpComplex_t> temp_1(N);
    temp_1.setZero();

    temp_1 = temp_matrix.colPivHouseholderQr().solve(temp_vec_1);

    mpComplex_t amplitude_1 {"0.0", "0.0"};

    for(int i = 0; i < N; ++i){
        amplitude_1 += abs(temp_1[i]);
    }
    amplitude_1 /= N;

    VectorX<mpComplex_t> stochastic_1(N);
    stochastic_1.setZero();

    for(int i = 0; i < N; ++i){
        stochastic_1[i] = amplitude_1 + (mpFloat_t(rand()) / mpFloat_t(RAND_MAX) +
                                         mpFloat_t("0.5")) * alpha;
    }

    VectorX<mpComplex_t> res_1(N);
    res_1.setZero();

    res_1 = temp_1 + stochastic_1;

    // II

    VectorX<mpComplex_t> temp_vec_2(N);
    temp_vec_2.setZero();

    for(int j = 0; j < N; ++j){
        temp_vec_2[j] = mpFloat_t("2.0") * PI * I *
                        (uInit[1] + rightHand1(uInit) * (l[j] - l_0));
    }

    VectorX<mpComplex_t> temp_2(N);
    temp_2.setZero();

    temp_2 = temp_matrix.colPivHouseholderQr().solve(temp_vec_2);

    mpComplex_t amplitude_2 {"0.0", "0.0"};

    for(int i = 0; i < N; ++i){
        amplitude_2 += abs(temp_2[i]);
    }
    amplitude_2 /= N;

    VectorX<mpComplex_t> stochastic_2(N);
    stochastic_2.setZero();

    for(int i = 0; i < N; ++i){
        stochastic_2[i] = amplitude_2 + (mpFloat_t(rand()) / mpFloat_t(RAND_MAX) +
                                         mpFloat_t("0.5")) * alpha;
    }

    VectorX<mpComplex_t> res_2(N);
    res_2.setZero();

    res_2 = temp_2 + stochastic_2;

    for(int j = 0; j < N; ++j){
        res[j] = res_1[j];
        res[j + N] = res_2[j];
    }

    this->uWaveInit = res;
}

void ivpCauchySolver::reCalc()
{
    VectorX<mpComplex_t> sol(2 * N);
    sol.setZero();
    VectorX<mpComplex_t> u_1(2 * N);
    u_1.setZero();
    VectorX<mpComplex_t> u_2(2 * N);
    u_2.setZero();

    u_1[0] = uInit[0];
    u_2[0] = uInit[1];

    u[0] = uInit[0];
    u[N] = uInit[1];

    mpComplex_t temp_1 {"0.0", "0.0"};
    mpComplex_t temp_2 {"0.0", "0.0"};

    for(int j = 1; j < N; ++j){
        temp_1 = mpComplex_t {"0.0", "0.0"};
        temp_2 = mpComplex_t {"0.0", "0.0"};

        for(int i = 0; i < N; ++i){
            temp_1 += uWave[i] * h_z[i] / (z[i] - l[j]);
            temp_2 += uWave[i + N] * h_z[i] / (z[i] - l[j]);
        }
        u_1[j] = real(temp_1 / (mpComplex_t("2.0") * PI * I));
        u[j] = u_1[j];
        u_2[j] = real(temp_2 / (mpComplex_t("2.0") * PI * I));
        u[j + N] = u_2[j];
    }
}

mpComplex_t ivpCauchySolver::rightHand0(VectorX<mpComplex_t> x)
{
    mpComplex_t temp {"0.0", "0.0"};
    temp = mpFloat_t("-1.0") * (pLambda0 * cos(x[1]) * (x[0] * x[0] - pA * pA) *
                                (x[0] * x[0] - pA * pA) / (x[0] * x[0] + pA * pA));
    return (temp / sqrt(mpFloat_t("1.0") + temp * temp));
}

mpComplex_t ivpCauchySolver::rightHand1(VectorX<mpComplex_t> x)
{
    mpComplex_t temp {"0.0", "0.0"};
    temp = mpFloat_t("-1.0") * (pLambda0 *
                                cos(x[1]) *
                                (x[0] * x[0] - pA * pA) *
                                (x[0] * x[0] - pA * pA) /
                                (x[0] * x[0] + pA * pA));
    return (mpFloat_t("1.0") / sqrt(mpFloat_t("1.0") + temp * temp));
}

mpComplex_t ivpCauchySolver::rightHandDerivative00(VectorX<mpComplex_t> x)
{
    mpComplex_t temp {"0.0", "0.0"};
    mpComplex_t tempU {"0.0", "0.0"};
    temp = mpFloat_t("-1.0") * (pLambda0 *
                                cos(x[1]) *
                                (x[0] * x[0] - pA * pA) *
                                (x[0] * x[0] - pA * pA) /
                                (x[0] * x[0] + pA * pA));
    tempU = mpFloat_t("-2.0") * (pLambda0 *
                                 cos(x[1]) *
                                 x[0] *
                                 (x[0] * x[0] - pA * pA) *
                                 (x[0] * x[0] + mpFloat_t("3.0") * pA * pA) /
                                 (x[0] * x[0] + pA * pA) /
                                 (x[0] * x[0] + pA * pA));
    return (tempU / ((mpFloat_t("1.0") + temp * temp) *
                      sqrt(mpFloat_t("1.0") + temp * temp)));
}

mpComplex_t ivpCauchySolver::rightHandDerivative01(VectorX<mpComplex_t> x)
{
    mpComplex_t temp {"0.0", "0.0"};
    mpComplex_t tempT {"0.0", "0.0"};
    temp = mpFloat_t("-1.0") * (pLambda0 *
                                cos(x[1]) *
                                (x[0] * x[0] - pA * pA) *
                                (x[0] * x[0] - pA * pA) /
                                (x[0] * x[0] + pA * pA));
    tempT = (pLambda0 *
             sin(x[1]) *
             (x[0] * x[0] - pA * pA) *
             (x[0] * x[0] - pA * pA) /
             (x[0] * x[0] + pA * pA));
    return (tempT / ((mpFloat_t("1.0") + temp * temp) *
                      sqrt(mpFloat_t("1.0") + temp * temp)));
}

mpComplex_t ivpCauchySolver::rightHandDerivative10(VectorX<mpComplex_t> x)
{
    mpComplex_t temp {"0.0", "0.0"};
    mpComplex_t tempU {"0.0", "0.0"};
    temp = mpFloat_t("-1.0") * (pLambda0 *
                                cos(x[1]) *
                                (x[0] * x[0] - pA * pA) *
                                (x[0] * x[0] - pA * pA) /
                                (x[0] * x[0] + pA * pA));
    tempU = mpFloat_t("-2.0") * (pLambda0 *
                                 cos(x[1]) *
                                 x[0] *
                                 (x[0] * x[0] - pA * pA) *
                                 (x[0] * x[0] + mpFloat_t("3.0") * pA * pA) /
                                 (x[0] * x[0] + pA * pA) /
                                 (x[0] * x[0] + pA * pA));
    return (mpFloat_t("-1.0") * temp * tempU / ((mpFloat_t("1.0") + temp * temp) *
                                                sqrt(mpFloat_t("1.0") + temp * temp)));
}

mpComplex_t ivpCauchySolver::rightHandDerivative11(VectorX<mpComplex_t> x)
{
    mpComplex_t temp {"0.0", "0.0"};
    mpComplex_t tempT {"0.0", "0.0"};
    temp = mpFloat_t("-1.0") * (pLambda0 *
                                cos(x[1]) *
                                (x[0] * x[0] - pA * pA) *
                                (x[0] * x[0] - pA * pA) /
                                (x[0] * x[0] + pA * pA));
    tempT = (pLambda0 *
             sin(x[1]) *
             (x[0] * x[0] - pA * pA) *
             (x[0] * x[0] - pA * pA) /
             (x[0] * x[0] + pA * pA));
    return (mpFloat_t("-1.0") * temp * tempT / ((mpFloat_t("1.0") + temp * temp) *
                                                sqrt(mpFloat_t("1.0") + temp * temp)));
}

VectorX<mpComplex_t> ivpCauchySolver::F(VectorX<mpComplex_t> x)
{
    VectorX<mpComplex_t> res (2 * N);
    res.setZero();

    mpComplex_t temp_1 {"0.0", "0.0"};
    mpComplex_t temp_2 {"0.0", "0.0"};

    for(int i = 0; i < N; ++i){
        temp_1 = temp_1 + x[i] * h_z[i] / (z[i] - l_0);
        temp_2 = temp_2 + x[i+N] * h_z[i] / (z[i] - l_0);
    }

    res[0] = mpFloat_t("2.0") * PI * I * uInit[0] - temp_1;
    res[N] = mpFloat_t("2.0") * PI * I * uInit[1] - temp_2;

    mpComplex_t temp_1_2 {"0.0", "0.0"};
    mpComplex_t temp_2_2 {"0.0", "0.0"};
    VectorX<mpComplex_t> tempV (2);
    tempV.setZero();

    for(int j = 1; j < N; ++j){
        temp_1 = mpComplex_t {"0.0", "0.0"};
        temp_1_2 = mpComplex_t {"0.0", "0.0"};
        temp_2 = mpComplex_t {"0.0", "0.0"};
        temp_2_2 = mpComplex_t {"0.0", "0.0"};
        for(int i = 0; i < N; ++i){
            temp_1 = temp_1 + x[i] * h_z[i] / (z[i] - l[j]);
            temp_1_2 = temp_1_2 + x[i] * h_z[i] / (z[i] - l[j]) / (z[i] - l[j]);
            temp_2 = temp_2 + x[i + N] * h_z[i] / (z[i] - l[j]);
            temp_2_2 = temp_2_2 + x[i + N] * h_z[i] / (z[i] - l[j]) / (z[i] - l[j]);
        }

        tempV[0] = temp_1 / (mpFloat_t("2.0") * PI * I);
        tempV[1] = temp_2 / (mpFloat_t("2.0") * PI * I);

        res[j] = mpFloat_t("2.0") * PI * I * rightHand0(tempV) - temp_1_2;
        res[j + N] = mpFloat_t("2.0") * PI * I * rightHand1(tempV) - temp_2_2;
    }
    return res;
}

MatrixX<mpComplex_t> ivpCauchySolver::JF(VectorX<mpComplex_t> x)
{
    MatrixX<mpComplex_t> res(2 * N, 2 * N);
    res.setZero();

    // 1

    for(int i = 0; i < N; ++i){
        res(0, i) = mpFloat_t("-1.0") * h_z[i] / (z[i] - l_0);
    }


    // 2 и 3

    mpComplex_t temp_1 {"0.0", "0.0"};
    mpComplex_t temp_2 {"0.0", "0.0"};
    VectorX<mpComplex_t> tempV(2);
    tempV.setZero();

    for(int j = 1; j < N; ++j){
        temp_1 = mpComplex_t{"0.0", "0.0"};
        temp_2 = mpComplex_t{"0.0", "0.0"};

        for(int k = 0; k < N; ++k){
            temp_1 = temp_1 + x[k] * h_z[k] / (z[k] - l[j]);
            temp_2 = temp_2 + x[k + N] * h_z[k] / (z[k] - l[j]);
        }

        tempV[0] = temp_1 / (mpFloat_t("2.0") * PI * I);
        tempV[1] = temp_2 / (mpFloat_t("2.0") * PI * I);

        for(int i = 0; i < N; ++i){
            res(j, i) = rightHandDerivative00(tempV) * h_z[i] / (z[i] - l[j]) -
                        h_z[i] / (z[i] - l[j]) / (z[i] - l[j]);
            res(j, i + N) = rightHandDerivative01(tempV) * h_z[i] / (z[i] - l[j]);
        }
    }

    // 4

    for(int i = 0; i < N; ++i){
        res(N, i + N) = mpFloat_t("-1.0") * h_z[i] / (z[i] - l_0);
    }

    // 5 и 6

    for(int j = 1; j < N; ++j){
        temp_1 = mpComplex_t{"0.0", "0.0"};
        temp_2 = mpComplex_t{"0.0", "0.0"};
        for(int k = 0; k < N; ++k){
            temp_1 = temp_1 + x[k] * h_z[k] / (z[k] - l[j]);
            temp_2 = temp_2 + x[k + N] * h_z[k] / (z[k] - l[j]);
        }

        tempV[0] = temp_1 / (mpFloat_t("2.0") * PI * I);
        tempV[1] = temp_2 / (mpFloat_t("2.0") * PI * I);

        for(int i = 0; i < N; ++i){
            res(j + N, i) = rightHandDerivative10(tempV) * h_z[i] / (z[i] - l[j]);
            res(j + N, i + N) = rightHandDerivative11(tempV) * h_z[i] / (z[i] - l[j]) -
                                h_z[i] / (z[i] - l[j]) / (z[i] - l[j]);
        }
    }
    return res;
}

void ivpCauchySolver::exactSolution(VectorX<mpComplex_t> x)
{
    VectorX<mpComplex_t> res(2 * N);
    res.setZero();

    for(int i = 0; i < 2 * N; ++i){
        if(i < N){
            res[i] = mpFloat_t("-2.0") * (pLambda0 * sin(x[i]) * pA * pA /
                                          (mpFloat_t("1.0") + sqrt(mpFloat_t("1.0") +
                                                                   mpFloat_t("4.0") * pLambda0 * pLambda0 *
                                                                       sin(x[i]) * sin(x[i]) * pA * pA)));
        }else{
            res[i] = x[i - N];
        }
    }

    uExact = res;
}

mpComplex_t ivpCauchySolver::relErr()
{
    mpComplex_t temp_1 {"0.0", "0.0"};
    mpComplex_t temp_2 {"0.0", "0.0"};
    mpComplex_t temp_d {"0.0", "0.0"};

    for(int j = 0; j < N - 1; ++j){
        temp_1 += (((uExact[j] - u[j]) * (uExact[j] - u[j]) +
                    (uExact[j + 1] - u[j + 1]) * (uExact[j + 1] - u[j + 1])) /
                   mpFloat_t("2.0"));
        temp_2 += (((uExact[j + N] - u[j + N]) * (uExact[j + N] - u[j + N]) +
                    (uExact[j + 1 + N] - u[j + 1 + N]) * (uExact[j + 1 + N] - u[j + 1 + N])) /
                   mpFloat_t("2.0"));
        temp_d += ((uExact[j] * uExact[j] +
                    uExact[j + 1] * uExact[j + 1] +
                    uExact[j + N] * uExact[j + N] +
                    uExact[j + 1 + N] * uExact[j + 1 + N]) /
                   mpFloat_t("2.0"));
    }

    return sqrt((temp_1 + temp_2) / temp_d);
}

void ivpCauchySolver::writeData(const int pr)
{
    ofs << "l;phi;z;h_z;uWaveInit;tWaveInit;uWave;tWave;u;t;uExact;tExact" << std::endl;
    for(int i = 0; i < N; ++i){
        // std::numeric_limits<mpComplex_t::value_type>::digits10
        ofs
            << mstr(l[i], pr) << ";"
            << mstr(phi[i], pr) << ";"
            << mstr(z[i], pr) << ";"
            << mstr(h_z[i], pr) << ";"
            << mstr(uWaveInit[i], pr) << ";"
            << mstr(uWaveInit[i + N], pr) << ";"
            << mstr(uWave[i], pr) << ";"
            << mstr(uWave[i + N], pr) << ";"
            << mstr(u[i], pr) << ";"
            << mstr(u[i + N], pr) << ";"
            << mstr(uExact[i], pr) << ";"
            << mstr(uExact[i + N], pr) << ";"
            << std::endl;
    }

    // *csvfile << "l" << "phi" << "z" << "h_z" << endrow;
    // *csvfile << "l" << endrow;
    // mpFloat_t temp_n;
    // string temp_s = "\0";
    // for(int i = 0; i < N; ++i){
        // temp_s = "\0";
        // temp_s += l[i].real().str();
        // temp_s += "+";
        // temp_s += l[i].imag().str();
        // temp_s += "i";
        // *csvfile << l[i].str() << phi[i].str() << z[i].str() << h_z[i].str() << endrow;
        // *csvfile << temp_s << endrow;
    // }
}

VectorX<mpComplex_t> ivpCauchySolver::testFunc(VectorX<mpComplex_t> x)
{
    Eigen::VectorX<mpComplex_t> res(2);
    res[0] = exp(x[0]) - x[1];
    res[1] = 1.0 - x[0] * x[0] - x[1] * x[1];
    return res;
}

MatrixX<mpComplex_t> ivpCauchySolver::testJFunc(VectorX<mpComplex_t> x)
{
    Eigen::MatrixX<mpComplex_t> res(2, 2);
    res(0, 0) = exp(x[0]);
    res(0, 1) = mpComplex_t(-1.0);
    res(1, 0) = mpComplex_t(-2.0) * x[0];
    res(1, 1) = mpComplex_t(-2.0) * x[1];
    return res;
}
