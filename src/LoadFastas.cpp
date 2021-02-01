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

#include "LoadFastas.h"
#include "Fasta.h"
#include "Main.h"

#include <QLabel>
#include <QLineEdit>
#include <QPushButton>
#include <QCheckBox>
#include <QVBoxLayout>

#include <hcsrc/FileReader.h>
#include <h3dsrc/Dialogue.h>

LoadFastas::LoadFastas(QWidget *parent) : QMainWindow(parent)
{
	QWidget *window = new QWidget();
	QVBoxLayout *box = new QVBoxLayout();
	window->setLayout(box);

	QHBoxLayout *hbox1 = new QHBoxLayout();
	QLabel *l = new QLabel("Fasta file: ", window);
	hbox1->addWidget(l);
	QLineEdit *e = new QLineEdit(window);
	_fastaLine = e;
	hbox1->addWidget(e);
	QPushButton *b = new QPushButton("Choose...", window);
	hbox1->addWidget(b);
	box->addLayout(hbox1);
	connect(b, &QPushButton::clicked, this, &LoadFastas::chooseFasta);

	QHBoxLayout *hbox2 = new QHBoxLayout();
	l = new QLabel("Focus on nucleotide range: ", window);
	hbox2->addWidget(l);
	e = new QLineEdit(window);
	_rangeLine = e;
	e->setPlaceholderText("x-y");
	hbox2->addWidget(e);
	box->addLayout(hbox2);

	QHBoxLayout *hbox4 = new QHBoxLayout();
	QCheckBox *c = new QCheckBox(window);
	_isProtein = c;
	l = new QLabel("Is protein sequence", window);
	hbox4->addWidget(c);
	hbox4->addWidget(l);
	box->addLayout(hbox4);

	QHBoxLayout *hbox3 = new QHBoxLayout();
	b = new QPushButton("Load", window);
	hbox3->addWidget(b);
	box->addLayout(hbox3);

	setCentralWidget(window);
	connect(b, &QPushButton::clicked, this, &LoadFastas::loadChosenFasta);
}

void LoadFastas::setProtein(bool p)
{
	_isProtein->setChecked(p);
}

void LoadFastas::loadChosenFasta()
{
	std::string file = _fastaLine->text().toStdString();
	std::string text = _rangeLine->text().toStdString();
	trim(text);
	std::vector<std::string> bits = split(text, '-');
	
	int start = -1;
	int end = -1;

	if (bits.size() == 2)
	{
		start = atoi(bits[0].c_str());
		end = atoi(bits[1].c_str());
	}

	loadFastas(file, start, end);
}

void LoadFastas::loadFastas(std::string filename, int start, int end)
{
	bool protein = (_isProtein->isChecked());
	loadSequence(filename, start, end, protein);
}

std::string makeSeq(std::vector<std::string> &lines, size_t *next)
{
	bool found_end = false;
	std::string seq = "";

	while (!found_end)
	{
		seq += lines[*next];
		
		while (seq.length() && (seq.back() < 'A' || seq.back() > 'Z'))
		{
			seq.pop_back();
		}

		(*next)++;
		
		if (*next >= lines.size() || lines[*next][0] == '>')
		{
			(*next)--;
			return seq;
		}
	}

	return "";
}

void LoadFastas::loadSequence(std::string filename, int start, int end,
                              bool isProtein)
{
	_fastaLine->setText("");

	if (!file_exists(filename))
	{
		std::cout << "no file like that" << std::endl;
		return;
	}

	std::string results = get_file_contents(filename);
	std::vector<std::string> lines = split(results, '\n');
	
	std::cout << "Focus: " << start << " " << end << std::endl;
	int count = 0;
	
	for (size_t i = 0; i < lines.size(); i++)
	{
		if (lines[i].length() <= 2 || lines[i][0] != '>')
		{
			continue;
		}

		std::string name = lines[i].substr(1);
		
		if (name.back() == '\n' || name.back() == '\r')
		{
			name.pop_back();
		}

		i++;
		std::string seq = makeSeq(lines, &i);
		
		if (!(start < 0 && end < 0))
		{
			if ((int)seq.length() < end)
			{
				if ((int)seq.length() < start)
				{
					continue;
				}

				seq = seq.substr(start);
			}
			else
			{
				seq = seq.substr(start, end - start);
			}
		}
		
		if (seq.length() == 0)
		{
			continue;
		}
		
		count++;
		Fasta *f = new Fasta(name);
		std::cout << count << ": Found " << name << std::endl;
		f->setSequence(seq, isProtein);
		
		_main->receiveSequence(f);
	}
	
	_main->makeSequenceMenu();
}

void LoadFastas::setMain(Main *m)
{
	_main = m;
}

void LoadFastas::chooseFasta()
{
	std::string filename = openDialogue(this, "Choose Fasta", 
	                                    "Fasta file (*.fasta)");

	_fastaLine->setText(QString::fromStdString(filename));
}
