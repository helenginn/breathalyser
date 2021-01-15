#include <cstdlib>
#include <iostream>
#include <QtWidgets/qapplication.h>
#include "Screen.h"
#include "Main.h"
#include "commit.h"

int main(int argc, char * argv[])
{
	std::cout << "Qt version: " << qVersion() << std::endl;

	QApplication app(argc, argv);
	setlocale(LC_NUMERIC, "C");
	srand(time(NULL));
	
	std::cout << "Check version: " << CHECK_VERSION_COMMIT_ID << std::endl;
	
	Main main(NULL);

	int status = app.exec();
	
	return status;
}
