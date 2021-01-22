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

#include "SlidingWindow.h"
#include "FastaMaster.h"
#include "Main.h"

#include <iostream>

#include <QVBoxLayout>
#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QCheckBox>
#include <QSlider>

SlidingWindow::SlidingWindow(QWidget *parent, Main *main) : QMainWindow(parent)
{
	QWidget *window = new QWidget();
	QVBoxLayout *box = new QVBoxLayout();
	window->setLayout(box);
	_main = main;

	FastaMaster *master = _main->fMaster();
	_scale = 100;
	int end = 20;
	
	if (master->fastaCount() < 1000)
	{
		_scale = 1;
		end = master->fastaCount();
	}

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
		QLabel *l = new QLabel("Window size in genomes:", window);
		hbox->addWidget(l);
		QSlider *s = new QSlider(Qt::Horizontal, window);
		int start = 20;
		s->setRange(0, end);
		s->setValue(start);

		connect(s, &QSlider::valueChanged,
		        this, &SlidingWindow::slid);

		_genomeCount = s;
		hbox->addWidget(s);
		QString startStr = QString::number(start * _scale);
		l = new QLabel(startStr, window);
		_genomeReport = l;
		hbox->addWidget(l);
		box->addLayout(hbox);
	}

	{
		QHBoxLayout *hbox = new QHBoxLayout();
		QCheckBox *c = new QCheckBox(window);
		c->setChecked(true);
		_overwrite = c;
		hbox->addWidget(c);
		QLabel *l = new QLabel("Overwrite existing frames", window);
		hbox->addWidget(l);
		box->addLayout(hbox);
	}
	
	{
		QHBoxLayout *hbox = new QHBoxLayout();
		QPushButton *b = new QPushButton("Run...", window);
		hbox->addWidget(b);
		box->addLayout(hbox);
		
		connect(b, &QPushButton::clicked, this, &SlidingWindow::run);
	}

	resize(400, 0);
	setCentralWidget(window);
}

void SlidingWindow::slid(int val)
{
	_genomeReport->setText(QString::number(val * _scale));
}

void SlidingWindow::run()
{
	int req = INT_MAX;
	unsigned char aa = '\0';
	size_t window_size = _genomeCount->value() * _scale;

	std::string str = _mutLine->text().toStdString();
	if (str.length())
	{
		char *ptr;
		req  = strtol(str.c_str(), &ptr, 0);
		
		if (ptr != NULL)
		{
			aa = *ptr;
		}
		
		std::cout << "Requiring residue " << req;
		
		if (aa > '\0')
		{
			std::cout << " to be " << aa << std::endl;
		}
		else
		{
			std::cout << " to be mutated" << std::endl;
		}
	}
	
	std::cout << "Window size: " << window_size << " genomes" << std::endl;
	
	bool overwrite = _overwrite->isChecked();
	
	std::cout << "Sliding window should " << (overwrite ? " " : "not ") << 
	"be overwriting existing frames." << std::endl;

	_main->slidingWindow(window_size, str, overwrite);
}

