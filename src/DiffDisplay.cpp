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

#include <QLabel>
#include <QPushButton>
#include <QSlider>
#include <QVBoxLayout>
#include <QTabWidget>

#include "Main.h"
#include "DiffDisplay.h"
#include "Difference.h"
#include "StructureView.h"

DiffDisplay::DiffDisplay(QWidget *parent, Difference *diff) : QWidget(parent)
{
	_diff = diff;

	QVBoxLayout *box = new QVBoxLayout();

	_label = new QLabel(this);
	_label->setScaledContents(true);
	_label->setSizePolicy(QSizePolicy::Preferred, QSizePolicy::Preferred);
	box->addWidget(_label);
	
	QHBoxLayout *hbox = new QHBoxLayout();
	QLabel *thLabel = new QLabel("Threshold:", this);
	_threshold = new QSlider(Qt::Horizontal, this);
	_threshold->setMinimum(0);
	_threshold->setMaximum(10000);
	_threshold->setTickInterval(1000);
	_threshold->setSingleStep(1);
	_threshold->setPageStep(100);
	_threshold->setTracking(true);
	hbox->addWidget(thLabel);
	hbox->addWidget(_threshold);
	
	_calculate = new QPushButton("Calculate", this);
	hbox->addWidget(_calculate);

	box->addLayout(hbox);
	
	setLayout(box);
}

void DiffDisplay::changeDifference(Difference *diff)
{
	if (diff != NULL)
	{
		_diff = diff;
	}
	
	if (_diff == NULL)
	{
		return;
	}

	_diff->setDisplay(this);
	int dim = _diff->atomCount() * 4;
	dim = std::min(dim, 800);
	QImage i = _diff->scaled(dim, dim, Qt::IgnoreAspectRatio);
	_label->setPixmap(QPixmap::fromImage(i));

	disconnect(_threshold, &QSlider::valueChanged, nullptr, nullptr);
	connect(_threshold, &QSlider::valueChanged,
	        _diff, &Difference::thresholdChanged);

	disconnect(_calculate, &QPushButton::clicked, nullptr, nullptr);
	connect(_calculate, &QPushButton::clicked,
	        _diff, &Difference::calculate);
	
	QTabWidget *tabs = _main->tabs();
	StructureView *view = _main->coupleView();

	connect(_calculate, &QPushButton::clicked,
	        tabs, [=]{ tabs->setCurrentWidget(view);});
}
