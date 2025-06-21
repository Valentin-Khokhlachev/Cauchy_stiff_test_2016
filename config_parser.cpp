#include "config_parser.h"

std::vector<std::string> parse_config_file(const std::string name, const char separator, const std::string key_to_find)
{
    using namespace std;

    vector<string> res;
    res.clear();

    string key = "";
    string tempString = "";
    char tempChar = '\0';

    bool flagFirstWord = true;
    bool flagKeyFound = false;
    bool flagKeyLineRead = false;
    bool flagKey = false;
    unsigned int counter = 0;

    ifstream ifs;
    ifs.open(name);

    if (ifs.is_open()) {
        while (ifs.get(tempChar) && !flagKey) {
            if (tempChar == '\n' || tempChar == '\0') {
                if (DEBUG_MOD) {
                    qDebug() << tempString.c_str();
                    qDebug() << "===";
                    qDebug() << key.c_str();
                    qDebug() << "======";
                }

                if (flagKeyFound) {
                    flagKeyLineRead = true;
                } else {
                    counter = 0;
                    tempString = "";
                    flagFirstWord = true;
                }
            } else {
                if (flagFirstWord && (tempChar == separator)) {
                    flagFirstWord = false;

                    key = tempString;
                    if (key == key_to_find) {
                        flagKeyFound = true;
                    }
                }

                if (DEBUG_MOD) {
                    qDebug() << ">>>" << tempChar << "<<<";
                }

                tempString += tempChar;
            }

            flagKey = flagKeyFound * flagKeyLineRead;

            if (tempChar == separator) {
                ++counter;
            }
        }

        if (flagKeyFound) {
            res.resize(counter);
            for (int i = 0; i < counter; ++i) {
                res[i] = "";
            }

            int res_counter = 0;
            int iterator = 0;
            flagFirstWord = true;
            tempChar = '\0';

            tempChar = tempString[0];
            while (tempChar != '\0') {
                tempChar = tempString[iterator];
                if (DEBUG_MOD) {
                    qDebug() << ">>" << tempChar << "<<";
                }

                if (tempChar && (tempChar != separator)) {
                    if (!flagFirstWord) {
                        if (DEBUG_MOD) {
                            qDebug() << ">" << tempChar << "<";
                        }
                        res[res_counter] += tempChar;
                    }
                } else {
                    if (flagFirstWord) {
                        flagFirstWord = false;
                    } else {
                        ++res_counter;
                    }
                }

                ++iterator;
            }
        } else {
            qDebug() << "There is no such key in \"" << name.c_str() << "\" !";
        }
        ifs.close();
    } else {
        qDebug() << "Can't open file!";
    }

    return res;
}

void myMessageHandler(QtMsgType type, const QMessageLogContext &, const QString & msg)
{
    QString txt;
    switch (type) {
    case QtDebugMsg:
        txt = QString("Debug: %1").arg(msg);
        break;
    case QtWarningMsg:
        txt = QString("Warning: %1").arg(msg);
        break;
    case QtCriticalMsg:
        txt = QString("Critical: %1").arg(msg);
        break;
    case QtFatalMsg:
        txt = QString("Fatal: %1").arg(msg);
        abort();
    }
    QFile outFile("../../../vpa_stiff_test_2016_arg_arc/results/log.txt");
    outFile.open(QIODevice::WriteOnly | QIODevice::Append);
    QTextStream ts(&outFile);
    ts << txt << Qt::endl;
}
