#include <QCoreApplication>
#include <QtCore>
#include <QVector>
#include <QDebug>
#include <QString>
#include <QtDebug>
#include <QFile>
#include <QTextStream>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <string>
#include <cassert>
#include <cstdlib>
#include <filesystem>
#include <memory>
#include <functional>

#include <eigen3/Eigen/Eigen>

#include "mpTypes.h"
#include "myTimer.h"
#include "csvstream.h"
#include "ivp_cauchy_solver.h"
#include "newton_iterator.h"
#include "config_parser.h"



int main(int argc, char *argv[])
{
    QCoreApplication app(argc, argv);

    QFile logFile;
    logFile.setFileName("../../../vpa_stiff_test_2016_arg_arc/results/log.txt");
    logFile.resize(0);
    qInstallMessageHandler(myMessageHandler);

    Timer::time_var start = time_now();

    Cauchy::ivpCauchySolver slv;
    slv.solve();
    slv.writeData(8);

    Timer::time_var finish = time_now();

    qDebug() << "---------";
    qDebug() << "Time taken" << time_taken(start, finish);
    qDebug() << "---------";

    // Сколько знаков выводить у чисел с расширенной разрядностью
    // int precForPrint = std::numeric_limits<mpComplex_t::value_type>::digits10;

    return 0;
}


