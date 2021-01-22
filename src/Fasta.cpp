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

#include "Fasta.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <hcsrc/Blast.h>
#include <hcsrc/FileReader.h>

Fasta::Fasta(std::string name)
{
	_problematic = false;
	_compared = false;
	_name = name;
	_orf = -1;
	_stop = -1;
	_offset = -1;
}

void Fasta::setSequence(std::string seq, bool protein)
{
	if (protein)
	{
		_result = seq;
		return;
	}

	_seq = seq;
}

void Fasta::description()
{
	std::cout << "Sequence of " << _seq.size() << " nucleotides." << std::endl;
}

bool Fasta::findNextORF()
{
	bool success = true;
	if (true || _orf < 0)
	{
		success = findNextATG();
	}
	
	if (!success)
	{
		return false;
	}

	success = findNextStopCodon();
	
	if (!success)
	{
//		return findNextATG();
	}
	
	return success;
}

bool Fasta::findNextATG()
{
	_orf++;
	bool success = false;

	while (_orf < (int)_seq.length() - 3)
	{
		if (_seq.substr(_orf, 3) == "ATG")
		{
			_stop = _seq.length() - 1;
			success = true;
			break;
		}
		
		_orf++;
	}
	
	if (success)
	{
//		std::cout << "ATG at index: "  << _orf << std::endl;
	}
	
	return (success);
}

bool Fasta::findNextStopCodon()
{
	_stop--;
	int remainder = ((_stop - _orf) % 3);
	_stop -= remainder;
	bool success = false;

	while (_stop > _orf + 3)
	{
		std::string tmp = _seq.substr(_stop, 3);
		if (tmp == "TAA" || tmp == "TAG" || tmp == "TGA")
		{
			success = true;
			break;
		}
		
		_stop -= 3;
	}
	
	return success;
}

std::string aaLetter(std::string str)
{
	if (str == "TTT" || str == "TTC")
	{
		return "F";
	}
	else if (str == "TTA" || str == "TTG" || str == "CTT" ||
	         str == "CTC" || str == "CTA" || str == "CTG")
	{
		return "L";
	}
	else if (str == "ATT" || str == "ATC" || str == "ATA")
	{
		return "I";
	}
	else if (str == "ATG")
	{
		return "M";
	}
	else if (str == "GTT" || str == "GTC" || str == "GTA"
	         || str == "GTG")
	{
		return "V";
	}
	else if (str == "TCT" || str == "TCC" || str == "TCA"
	         || str == "TCG" || str == "AGT" || str == "AGC")
	{
		return "S";
	}
	else if (str == "CCT" || str == "CCC" || str == "CCA"
	         || str == "CCG")
	{
		return "P";
	}
	else if (str == "ACT" || str == "ACC" || str == "ACA"
	         || str == "ACG")
	{
		return "T";
	}
	else if (str == "GCT" || str == "GCC" || str == "GCA"
	         || str == "GCG")
	{
		return "A";
	}
	else if (str == "TAT" || str == "TAC")
	{
		return "Y";
	}
	else if (str == "CAT" || str == "CAC")
	{
		return "H";
	}
	else if (str == "CAA" || str == "CAG")
	{
		return "Q";
	}
	else if (str == "AAT" || str == "AAC")
	{
		return "N";
	}
	else if (str == "AAA" || str == "AAG")
	{
		return "K";
	}
	else if (str == "GAT" || str == "GAC")
	{
		return "D";
	}
	else if (str == "GAA" || str == "GAG")
	{
		return "E";
	}
	else if (str == "TGT" || str == "TGC")
	{
		return "C";
	}
	else if (str == "TGG")
	{
		return "W";
	}
	else if (str == "AGA" || str == "AGG" || str == "CGT"
	         || str == "CGC" || str == "CGA" || str == "CGG")
	{
		return "R";
	}
	else if (str == "GGA" || str == "GGG" || str == "GGT"
	         || str == "GGC")
	{
		return "G";
	}

	return " ";
}

std::string Fasta::generateSequence()
{
	std::string aaSeq = "";
	for (int i = _orf; i < _stop + 3; i += 3)
	{
		std::string tmp = _seq.substr(i, 3);
		
		std::string aa = aaLetter(tmp);
		aaSeq += aa;
	}
	
	return aaSeq;
}

bool Fasta::roughlyAlign(std::string mine, std::string ref, int minRes)
{
	if (mine == "" && hasResult())
	{
		mine = result();
	}
	else if (mine == "")
	{
		return false;
	}
	
	int size = 10;

	for (size_t i = 0; i < ref.length(); i += size / 2)
	{
		std::string probe = ref.substr(i, size);

		if (probe.length() < size)
		{
			break;
		}

		if (probe.find(" ") != std::string::npos)
		{
			continue;
		}

		/* position of probe in my sequence */
		size_t loc = mine.find(probe);

		if (loc != std::string::npos)
		{
			/* how much further forward we are in relation to ref */
			int diff = (int)loc - (int)i;
			
			/* which is also the beginning...*/
			int beginning = diff;
			
			/* but we have minRes residues before then, which we may want */
			beginning -= minRes;
			
			/* but this might be too much */
			if (beginning < 0)
			{
				minRes -= beginning;
				beginning = 0;
			}
			
			/* we do want the full length of the reference */
			int end = ref.length();
			
			/* but this might be too much */
			
			if ((int)end - (int)beginning > mine.length())
			{
				end = mine.length() - beginning;
			}

			_result = mine.substr(beginning, end);
			_offset = minRes;
			_seq.clear();

			std::cout << "Found " << name() << ", offset " << _offset << std::endl;
			return true;
		}
	}
	
	return false;
}

std::string Fasta::roughCompare(std::string ref, int minRes)
{
	if (hasResult())
	{
		if (roughlyAlign(_result, ref, minRes))
		{
			return _result;
		}

		return _result;
	}

	while (findNextORF())
	{
		std::string seq = generateSequence();
		size_t size = 10;
		
		if (roughlyAlign(seq, ref, minRes))
		{
			return seq;
		}
	}

	return "";
}

int Fasta::gapCount()
{
	int count = 0;
	char ch = ' ';
	for (size_t i = 0; (i = _result.find(ch, i)) != std::string::npos; i++) 
	{
		count++;
	}

	return count;
}

double Fasta::compareWithFasta(Fasta *f)
{
	if (_result.length() != f->_result.length())
	{
		return 0;
	}

	_mutations.clear();

	bool last = true;
	int best_run = 0;
	int run = 1;
	int muts = 0;

	for (size_t i = 0; i < _result.length(); i++)
	{
		char a = _result[i];
		char b = f->_result[i];
		
		if (a == ' ' || b == ' ' || a == b)
		{
			last = false;
			if (run > best_run)
			{
				best_run = run;
			}

			run = 1;

			continue;
		}
		
		if (last)
		{
			run++;
		}

		last = true;

		addMutation(b, i, a);

		muts++;
	}

	double distance = muts;
	distance /= _result.length() / 200;

	double score = exp(-(distance * distance));

	return score;
}

void Fasta::addMutation(char fromWhat, int mut, char towhat)
{
	std::ostringstream ss;
	ss << fromWhat << mut + _offset << towhat;
	std::string str = ss.str();

	_mutations.push_back(str);
}

void Fasta::carefulCompareWithFasta(Fasta *f)
{
	if (!hasResult() && !f->hasResult())
	{
		return;
	}

	carefulCompareWithString(f->result());
}

void Fasta::carefulCompareWithString(std::string seq2)
{
	std::string seq1 = result();
	
	_mutations.clear();
	
	if (seq1.length() == 0 || seq2.length() == 0)
	{
		return;
	}

	Alignment ala, alb;
	int muts, dels;
	srand(1);
	compare_sequences_and_alignments(seq1, seq2, &muts, &dels, ala, alb);
	tidy_alignments(ala, alb);

	std::ostringstream ssleft, ssalign, ssright;
	std::vector<int> indices;
	print_alignments(ala, alb, ssleft, ssalign, ssright, indices);
	
	std::string left =   ssleft.str();
	std::string align = ssalign.str();
	std::string right = ssright.str();
	
	/*
	std::cout << name() << std::endl;
	std::cout << left << std::endl;
	std::cout << align << std::endl;
	std::cout << right << std::endl;
	std::cout << std::endl;
	*/
	
	int running = 0;
	int offset = 0;
	int plus = 0;

	for (size_t i = 0; i < align.size(); i++)
	{
		if (align[i] == '.' || align[i] == '*')
		{
			offset += running;
			running = 0;
		}

		if (align[i] != '+' && !(left[i] == ' ' && right[i] == ' '))
		{
			plus = 0;
		}

		offset = 0;

		if (align[i] == '.')
		{
			continue;
		}

		if (align[i] == '*')
		{
			addMutation(right[i], indices[i] + offset, left[i]);
			continue;
		}
		
		if (align[i] == '+')
		{
			addMutation(left[i], indices[i] + offset, '+');
			running--;
			plus++;
		}
		
		if (align[i] == '-')
		{
			addMutation(right[i], indices[i] + offset, '-');
			running++;
		}
		
		if (plus > 6)
		{
			_problematic = true;
		}
	}
	
	_compared = true;
}

std::string Fasta::mutationSummary()
{
	std::ostringstream ss;
	
	for (size_t i = 0; i < mutationCount(); i++)
	{
		ss << mutation(i) << " ";
	}

	std::string str = ss.str();
	
	if (str.length() > 0)
	{
		str.pop_back();
	}

	return str;
}

void Fasta::loadMutations(std::string muts)
{
	_mutations.clear();
	
	std::vector<std::string> components = split(muts, ' ');
	
	for (size_t i = 0; i < components.size(); i++)
	{
		if (components[i].length() == 0)
		{
			continue;
		}
		
		_mutations.push_back(components[i]);
	}
	
	_compared = true;
}
