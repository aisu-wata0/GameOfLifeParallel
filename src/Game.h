#ifndef GAME_H
#define GAME_H

#include<iostream>
#include<unistd.h>
#include<stdlib.h>
#include<vector>
#include<string>
#include<sstream>

namespace G {

using namespace std;

#define LIVE 'O'
#define DEAD '.'

class Game
{
private:
	long size_;
	long sizeMem_;
	vector<char> grid_;
	vector<char> temp_;
	
public:
	
	Game(long size)
	{
		size_ = size;
		sizeMem_ = size_ + 2;
		
		grid_.reserve(sizeMem_*sizeMem_);
		temp_.reserve(sizeMem_*sizeMem_);
		
		for(long i = 0; i < sizeMem_; ++i){
			for(long j = 0; j < sizeMem_; ++j)
			{
				grid_.push_back(DEAD);
				temp_.push_back(DEAD);
			}
		}
	}
	
	long size() const {
		return size_;
	}
	
	char& cell(long i, long j){
		return grid_[(i+1)*sizeMem_ + (j+1)];
	}
	const char& cell(long i, long j) const {
		return grid_[(i+1)*sizeMem_ + (j+1)];
	}
	
	char& cellTmp(long i, long j){
		return temp_[(i+1)*sizeMem_ + (j+1)];
	}
	
	void print() const
	{
		for(long i = 0; i < size_; ++i){
			for(long j = 0; j < size_; ++j)
			{
				cout << cell(i,j) << " ";
			}
			cout << "\n";
		}
		cout << flush;
	}

	void nextGen()
	{
		
		// blocking
		/*
		long blockSZ = 80;

		for( long ib = 0; ib < size_; ib += blockSZ ){
		long max_i = std::min(ib + blockSZ, size_);
		for( long jb = 0; jb < size_; jb += blockSZ ){
		long max_j = std::min(jb + blockSZ, size_);
			for( long i = ib; i < max_i; ++i )
			for( long j = jb; j < max_j; ++j )
			{
			do_something(i,j);
			}
		}
		}
		*/
		
		for(long i = 0; i < size_; ++i)
		{
			for(long j = 0; j < size_; ++j)
			{
				if(cell(i,j)==LIVE)
				{
					if(neighbours(i,j)<2)
					{
						cellTmp(i,j) = DEAD;
					}
					else if(neighbours(i,j)>3)
					{
						cellTmp(i,j) = DEAD;
					}
					else
					{
						cellTmp(i,j) = LIVE;
					}
				}
				else
				{
					if(neighbours(i,j)==3)
					{
						cellTmp(i,j) = LIVE;
					}
					else
					{
						cellTmp(i,j) = DEAD;
					}
				}
			}
		}
		
		grid_.swap(temp_);
	}

	int neighbours(long x, long y) const
	{
		int count = 0;
		
		if(cell(x-1,y-1) == LIVE){
			++count;
		}
		if(cell(x-1,y) == LIVE){
			++count;
		}
		if(cell(x-1,y+1) == LIVE){
			++count;
		}
		if(cell(x,y-1) == LIVE){
			++count;
		}
		if(cell(x,y+1) == LIVE){
			++count;
		}
		if(cell(x+1,y-1) == LIVE){
			++count;
		}
		if(cell(x+1,y) == LIVE){
			++count;
		}
		if(cell(x+1,y+1) == LIVE){
			++count;
		}

		return(count);
	}
	
	void init(string s)
	{
		if(s == "Glider"){
			cell(0,1) = LIVE;
			cell(1,2) = LIVE;
			cell(2,0) = LIVE;
			cell(2,1) = LIVE;
			cell(2,2) = LIVE;
		}

		if(s == "Pulsar"){
			cell(2,4) = LIVE;
			cell(2,5) = LIVE;
			cell(2,6) = LIVE;
			cell(2,10) = LIVE;
			cell(2,11) = LIVE;
			cell(2,12) = LIVE;
			cell(4,2) = LIVE;
			cell(4,7) = LIVE;
			cell(4,9) = LIVE;
			cell(4,14) = LIVE;
			cell(5,2) = LIVE;
			cell(5,7) = LIVE;
			cell(5,9) = LIVE;
			cell(5,14) = LIVE;
			cell(6,2) = LIVE;
			cell(6,7) = LIVE;
			cell(6,9) = LIVE;
			cell(6,14) = LIVE;
			cell(7,4) = LIVE;
			cell(7,5) = LIVE;
			cell(7,6) = LIVE;
			cell(7,10) = LIVE;
			cell(7,11) = LIVE;
			cell(7,12) = LIVE;
			cell(9,4) = LIVE;
			cell(9,5) = LIVE;
			cell(9,6) = LIVE;
			cell(9,10) = LIVE;
			cell(9,11) = LIVE;
			cell(9,12) = LIVE;
			cell(10,2) = LIVE;
			cell(10,7) = LIVE;
			cell(10,9) = LIVE;
			cell(10,14) = LIVE;
			cell(11,2) = LIVE;
			cell(11,7) = LIVE;
			cell(11,9) = LIVE;
			cell(11,14) = LIVE;
			cell(12,2) = LIVE;
			cell(12,7) = LIVE;
			cell(12,9) = LIVE;
			cell(12,14) = LIVE;
			cell(14,4) = LIVE;
			cell(14,5) = LIVE;
			cell(14,6) = LIVE;
			cell(14,10) = LIVE;
			cell(14,11) = LIVE;
			cell(14,12) = LIVE;
		}

		if(s == "Glider gun"){
			cell(5,1) = LIVE;
			cell(6,1) = LIVE;
			cell(5,2) = LIVE;
			cell(6,2) = LIVE;

			cell(5,11) = LIVE;
			cell(6,11) = LIVE;
			cell(7,11) = LIVE;
			cell(4,12) = LIVE;
			cell(8,12) = LIVE;
			cell(3,13) = LIVE;
			cell(9,13) = LIVE;
			cell(3,14) = LIVE;
			cell(9,14) = LIVE;
			cell(6,15) = LIVE;
			cell(4,16) = LIVE;
			cell(8,16) = LIVE;
			cell(5,17) = LIVE;
			cell(6,17) = LIVE;
			cell(7,17) = LIVE;
			cell(6,18) = LIVE;

			cell(3,21) = LIVE;
			cell(4,21) = LIVE;
			cell(5,21) = LIVE;
			cell(3,22) = LIVE;
			cell(4,22) = LIVE;
			cell(5,22) = LIVE;
			cell(2,23) = LIVE;
			cell(6,23) = LIVE;
			cell(1,25) = LIVE;
			cell(2,25) = LIVE;
			cell(6,25) = LIVE;
			cell(7,25) = LIVE;

			cell(3,35) = LIVE;
			cell(4,35) = LIVE;
			cell(3,36) = LIVE;
			cell(4,36) = LIVE;
		}
	}
};

}
#endif // GAME_H
