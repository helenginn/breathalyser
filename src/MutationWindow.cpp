// breathalyser
// Copyright (C) 2019 Helen Ginn
// 
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
// 
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
// 
// Please email: vagabond @ hginn.co.uk for more details.

#include "MutationWindow.h"
#include "FastaMaster.h"
#include "Main.h"

#include <QVBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>

MutationWindow::MutationWindow(QWidget *parent) : QMainWindow(parent)
{
	QWidget *window = new QWidget();
	QVBoxLayout *box = new QVBoxLayout();
	window->setLayout(box);

	{
		QHBoxLayout *hbox = new QHBoxLayout();
		QLabel *l = new QLabel("Require mutation:", window);
		hbox->addWidget(l);
		QLineEdit *e = new QLineEdit(window);
		_mutLine = e;
		e->setPlaceholderText("e.g. 501 or 501Y");
		hbox->addWidget(e);
		box->addLayout(hbox);
	}

	{
		QHBoxLayout *hbox = new QHBoxLayout();
		QPushButton *b = new QPushButton("Set", window);
		hbox->addWidget(b);
		box->addLayout(hbox);
		
		connect(b, &QPushButton::clicked, this, &MutationWindow::run);
	}

	resize(400, 0);
	setCentralWidget(window);
}

void MutationWindow::run()
{
	std::string str = _mutLine->text().toStdString();
	FastaMaster *master = _main->fMaster();
	
	master->requireMutation(str);
	hide();
}

