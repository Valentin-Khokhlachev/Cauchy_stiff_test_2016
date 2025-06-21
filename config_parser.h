#ifndef CONFIG_PARSER_H
#define CONFIG_PARSER_H

#include <QCoreApplication>
#include <QtCore>
#include <QtDebug>
#include <QFile>
#include <QTextStream>
#include <QDebug>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "debug_mod.h"

std::vector<std::string> parse_config_file(const std::string name, const char separator, std::string key_to_find);
void myMessageHandler(QtMsgType type, const QMessageLogContext &, const QString & msg);

#endif // CONFIG_PARSER_H
