#include <utils/utils.h>


std::vector<Logger::ConsoleText> Logger::consoleOutput;
int Logger::maxConsoleLineCount = 15;

std::string Logger::ms_strLogPath = std::string(SCP_DATA_FOLDER "/out");
std::string Logger::ms_strPrintFileName =
    Logger::ms_strLogPath + std::string("/print.txt");
std::string Logger::ms_strLogFileName =
    Logger::ms_strLogPath + std::string("/log.txt");
std::string Logger::ms_strConsoleFileName =
    Logger::ms_strLogPath + std::string("/console.txt");

void getCharSeparatedStringList(const char *cString,
                                std::vector<std::string> &lines,
                                char separator) {
    lines.clear();
    std::string line;
    for (const char *c = cString; *c; c++) {
        if (*c && *c != separator)
            line.append(1, *c);
        else {
            lines.push_back(line);
            line.clear();
        }
    }
    if (line.size() > 0) lines.push_back(line);
}

void Logger::print(PRINT_COLOR color, char *pBuffer) {
    createPath(ms_strLogPath);

    // Write to file
    static FILE *fp = fopen(ms_strPrintFileName.c_str(), "wt");
    fprintf(fp, "%s", pBuffer);
    fflush(fp);

    // Write to terminal
    if (color == Logger::RED)
        printf("%s%s%s", ANSI_COLOR_RED, pBuffer, ANSI_COLOR_RESET);
    else if (color == Logger::GREEN)
        printf("%s%s%s", ANSI_COLOR_GREEN, pBuffer, ANSI_COLOR_RESET);
    else if (color == Logger::YELLOW)
        printf("%s%s%s", ANSI_COLOR_YELLOW, pBuffer, ANSI_COLOR_RESET);
    else if (color == Logger::BLUE)
        printf("%s%s%s", ANSI_COLOR_BLUE, pBuffer, ANSI_COLOR_RESET);
    else if (color == Logger::MAGENTA)
        printf("%s%s%s", ANSI_COLOR_MAGENTA, pBuffer, ANSI_COLOR_RESET);
    else if (color == Logger::CYAN)
        printf("%s%s%s", ANSI_COLOR_CYAN, pBuffer, ANSI_COLOR_RESET);
    else if (color == Logger::RESET)
        printf("%s%s%s", ANSI_COLOR_RESET, pBuffer, ANSI_COLOR_RESET);

    RELEASE_STRING_FROM_ARGUMENT_LIST(pBuffer);
}

void Logger::print(const char *fmt, ...) {
    char *pBuffer = nullptr;
    GET_STRING_FROM_ARGUMENT_LIST(fmt, pBuffer);
    print(RESET, pBuffer);
}

void Logger::print(Logger::PRINT_COLOR color, const char *fmt, ...) {
    char *pBuffer = nullptr;
    GET_STRING_FROM_ARGUMENT_LIST(fmt, pBuffer);
    print(color, pBuffer);
}

void Logger::logPrint(const char *fmt, ...) {
    char *pBuffer = nullptr;
    GET_STRING_FROM_ARGUMENT_LIST(fmt, pBuffer);

    createPath(ms_strLogPath);

    static FILE *fp = fopen(ms_strLogFileName.c_str(), "wt");
    fprintf(fp, "%s", pBuffer);
    fflush(fp);

    RELEASE_STRING_FROM_ARGUMENT_LIST(pBuffer);
}
void Logger::consolePrint(const Eigen::Vector3d &color, char *pBuffer) {
    createPath(ms_strLogPath);

    //Write to file
    static FILE *fp = fopen(ms_strConsoleFileName.c_str(), "wt");
    fprintf(fp, "%s", pBuffer);
    fflush(fp);

    //Store for console print
    std::vector<std::string> newLines;
    getCharSeparatedStringList(pBuffer, newLines);

    for (unsigned int i = 0; i < newLines.size(); i++)
        consoleOutput.push_back(ConsoleText{newLines[i], color});

    while ((int)consoleOutput.size() > maxConsoleLineCount)
        consoleOutput.erase(consoleOutput.begin());

    RELEASE_STRING_FROM_ARGUMENT_LIST(pBuffer);
}

void Logger::consolePrint(const char *fmt, ...) {
    char *pBuffer = nullptr;
    GET_STRING_FROM_ARGUMENT_LIST(fmt, pBuffer);
    consolePrint(Eigen::Vector3d::Ones(), pBuffer);
}

void Logger::consolePrint(const Eigen::Vector3d &color, const char *fmt, ...) {
    char *pBuffer = nullptr;
    GET_STRING_FROM_ARGUMENT_LIST(fmt, pBuffer);
    consolePrint(color, pBuffer);
}

FILE *Logger::getConsoleFilePointer() {
    static FILE *fp = fopen(ms_strConsoleFileName.c_str(), "wt");
    return fp;
}

bool Logger::createPath(const std::string &_strPath) {
    bool bSuccess = false;

#ifdef WIN32
    std::string stemp = std::string(_strPath.begin(), _strPath.end());
    LPCSTR sw = stemp.c_str();

    if (CreateDirectory(sw, nullptr)) {
        // Directory created
        bSuccess |= true;
    } else if (ERROR_ALREADY_EXISTS == GetLastError()) {
        // Directory already exists
        bSuccess |= true;
    }

#else
    bSuccess = (mkdir(_strPath.c_str(), 0777) != -1);
#endif
    return bSuccess;
}

Logger::Logger() {}

Logger::~Logger() {}
