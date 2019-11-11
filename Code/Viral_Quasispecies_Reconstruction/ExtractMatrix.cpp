#include <algorithm>
#include <random>
#include <iomanip> // std::setprecision
#include <vector> // for use of vector
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream> 
#include <numeric>
#include <time.h>

using namespace std;

// Code for parsing SAM file (line 59-339) adapted from PredictHaplo
// http://bmda.cs.unibas.ch/software.html

string binary(int number, stringstream& strs) 
{
	int remainder;
	if(number <= 1) 
	{
		strs << number;
	  	return strs.str() ;
	}
	remainder = number%2;
	binary(number >> 1, strs);    
	strs << remainder;
	
	return  strs.str();
}

vector<string> tokenize(const string& str,const string& delimiters)
{
	vector<string> tokens;	
	string::size_type lastPos = 0, pos = 0;  
  	int count = 0;
  	if(str.length()<1)  return tokens;
   
  	lastPos = str.find_first_not_of(delimiters, 0); // skip delimiters at beginning.     
  	if((str.substr(0, lastPos-pos).length()) > 0)
  	{	  	
  		count = str.substr(0, lastPos-pos).length();  	
  		for(int i=0; i < count; i++)  	
  	 		tokens.push_back("");
  		if(string::npos == lastPos)
	  		tokens.push_back("");
	}

  	pos = str.find_first_of(delimiters, lastPos); // find first \"non-delimiter\".
  	while (string::npos != pos || string::npos != lastPos)
 	 {  	      	    
     		tokens.push_back( str.substr(lastPos, pos - lastPos)); // found a token, add it to the vector.				
     		lastPos = str.find_first_not_of(delimiters, pos); // skip delimiters.  Note the \"not_of\"   	   	    
		
		if((string::npos != pos) && (str.substr(pos, lastPos-pos).length() > 1))  		
  		{
  			count = str.substr(pos, lastPos-pos).length();
  			for(int i=0; i < count-1; i++)
  	 			tokens.push_back("");
		}	
  		pos = str.find_first_of(delimiters, lastPos);
	}
	return tokens;
}

void parse_sam_line(const vector<string>& tokens, vector<int>& seq_b, double& a_score, char gap_quality, int& indels, int& al_start)
{  
	seq_b.clear();
	indels = 0;
	
  	bool foundAscore = false;
  	for(int j = 11; j< tokens.size();j++)
	{
  		vector<string> tokens_as = tokenize( tokens[j].c_str(),":");
    		if( tokens_as[0] == "AS")
		{
      			a_score = atof(tokens_as[2].c_str());
      			foundAscore = true;
      			break;
    		}
 	} 
	if (!foundAscore)
		cout << foundAscore << " no a_score "<< a_score << endl;
 
	al_start = atoi( tokens[3].c_str()) ;
 
	vector<int> isAlpha(tokens[5].size(),1);
  	for(int i =0; i< tokens[5].size();i++)
	{
    		if(!isalpha(tokens[5][i]))
      			isAlpha[i] = 0;
	}
  
  	vector<int> sub_length_vec;
  	vector<char> symbols;
  	int sub_length =0;
  	for(int i =0; i< tokens[5].size();i++)
	{
    		if(isAlpha[i] == 1)
		{
      			sub_length_vec.push_back( atoi(tokens[5].substr(i-sub_length,sub_length).c_str()));
      			symbols.push_back(tokens[5][i]);
      			sub_length =0;
    		}
    		if(isAlpha[i] == 0)
     			sub_length++;	
	}
  
  	int c =0;
  	for(int i =0; i< sub_length_vec.size();i++)
	{  
    		if(symbols[i] == 'S')
      			for(int j =0; j< sub_length_vec[i];j++)
				c++;
    		else if(symbols[i] == 'M')
			for(int j =0; j< sub_length_vec[i];j++)
			{
				int k = 0;
				if(tokens[9][c] == 'A' || tokens[9][c] == 'a')
			  		k=1;
				else if(tokens[9][c] == 'C' || tokens[9][c] == 'c')
			  		k=2;
				else if(tokens[9][c] == 'G' || tokens[9][c] == 'g')
			  		k=3;
				else if(tokens[9][c] == 'T' || tokens[9][c] == 't')
			  		k=4;
				seq_b.push_back(k);
				c++;
		      }
		else if(symbols[i] == 'I')
      			for(int j =0; j< sub_length_vec[i];j++)
			{	
				c++;
      			}
    		else if(symbols[i] == 'D')
      			for(int j =0; j< sub_length_vec[i];j++)
			{
				seq_b.push_back(0);
				indels++;	
      			}
	}  
}

int parseSAMpaired( string al, double min_qual, int min_length, int max_insertln, char gap_quality, double& mean_length, vector<vector<int> >& Read_matrix, int reconstruction_start, int reconstruction_end, int& total_count, int gene_length, string zonename, vector<vector<int> >& ReadSeq, vector<int>& StartSeq, vector<int>& EndSeq)
{
	string line;
  	vector<string> tokens, tokens_1, tokens_2;
	int pair_counter = 0, seq_counter =0, singleton_counter = 0;
  	vector<int> seq_b, seq_b_pairs;
  	double a_score,a_score_pairs;
	int indels, indels_pairs;
  	int al_start, al_start_pairs;
    	string sRC_1, sRC_2, pairs_1, pairs_2;
  	bool part_1 = false, is_pair = false;
  	string id, id_1;
  	int RC;
	vector<vector<int> > lowQSseq;
  
	int filtered_counter1 = 0, filtered_counter2 = 0;
	int mapped_counter = 0, unmapped_counter = 0; 
int filtered_cond2 = 0;

	ifstream inf6(al.c_str(),ios::in);
  	while( getline(inf6,line,'\n') )
	{
    		if(line[0] == '@')
      			continue;
	    	tokens = tokenize(line,"\t");
		if(tokens.size() <5)
		{
			cout << "Problem with sam file..." << endl;
			return 1;
		}
		id =  tokens[0];
	    	total_count++;
				
	    	RC =  atoi( tokens[1].c_str());
	    	stringstream strs;
	    	string sRC = binary(RC,  strs);
		int sz = sRC.size();
		if(sz > 8) // bit9-12 should be 0
	      	{ filtered_counter1++;	continue;}
	        if(sRC[sz-3] == '1' ||  sRC[sz-4]  == '1'  ) // bit 3-4 should be 0
	      	{ filtered_counter2++;	continue;}
	    	if(part_1 && id == id_1 && sz == 8 && sRC[0] == '1')
		{
	      		is_pair = true;
	      		part_1 = false;
	      		pairs_2 = line;
	      		sRC_2 = sRC;
	      		tokens_2 = tokens;
	      		pair_counter++;
	   	}
	    	else
		{
	     		part_1 = true;
	      		is_pair = false;
	      		pairs_1 = line;
	      		sRC_1 = sRC;
	      		id_1 = id;
	      		tokens_1 = tokens;
	      		singleton_counter++;
	   	}
		
	    	if(is_pair)
		{	
	      		parse_sam_line(tokens_1, seq_b,  a_score,  gap_quality, indels, al_start );
	       		parse_sam_line(tokens_2, seq_b_pairs,  a_score_pairs,  gap_quality, indels_pairs, al_start_pairs );
		
	      		int qual = atoi(tokens_1[4].c_str());
	      		int qual_pairs = atoi(tokens_2[4].c_str());

	      		if(// seq_b.size() >= min_length && seq_b_pairs.size() >= min_length 
//				&& seq_b.size() !=0 && seq_b_pairs.size() !=0 && 
				qual >= min_qual && qual_pairs >= min_qual)
			{	 
				mapped_counter++;
				int StartPos = al_start;
				vector<int> SEQ_combined = seq_b;
				bool is_gap = false;
				int Nlength = 0;
	
				if(al_start_pairs < al_start)
				{
		  			StartPos = al_start_pairs;	  
		  			is_gap = false;
		  			if(al_start-(al_start_pairs +   (int)seq_b_pairs.size())>0)
		    				is_gap = true;
	  				if(is_gap)
					{
						Nlength = al_start-(al_start_pairs +   (int)seq_b_pairs.size());
	    					vector<int> Ns(Nlength,0);
						seq_b_pairs.insert(seq_b_pairs.end(), Ns.begin(), Ns.end());
						seq_b_pairs.insert(seq_b_pairs.end(), seq_b.begin(), seq_b.end());
						SEQ_combined = seq_b_pairs; 
	  				}
					else
					{
						vector<int> first_part (seq_b_pairs.begin(),seq_b_pairs.begin() + (al_start - al_start_pairs));
						first_part.insert(first_part.end(), seq_b.begin(), seq_b.end());
						SEQ_combined = first_part;
	  				}
				}
				else
				{
					if(al_start_pairs > al_start)
					{
				    		is_gap = false;
				    		if(al_start_pairs - (al_start + (int)seq_b.size()) >0)
				      			is_gap = true;
	    					if(is_gap)
						{
	      						Nlength =al_start_pairs-(al_start + (int)seq_b.size());
	      						vector<int> Ns(Nlength,0);
	      						seq_b.insert(seq_b.end(), Ns.begin(), Ns.end());
	      						seq_b.insert(seq_b.end(), seq_b_pairs.begin(), seq_b_pairs.end());
	      						SEQ_combined = seq_b;
	   		 			}
						else
						{
	      						vector<int> first_part (seq_b.begin(),seq_b.begin() + (al_start_pairs - al_start));
	      						first_part.insert(first_part.end(), seq_b_pairs.begin(), seq_b_pairs.end());
	      						SEQ_combined = first_part;
	   					 }
	  				}
				}
				if(!is_gap || Nlength < max_insertln)  
				{filtered_cond2++;
					int EndPos = StartPos + SEQ_combined.size()-1; //range = reconstruction_end-reconstruction_start+1;
					if ( StartPos <= reconstruction_end && EndPos >= reconstruction_start)
					{
						vector<int> SEQ_range;
						if (StartPos < reconstruction_start)
						{
							if (EndPos <= reconstruction_end)
							{
								vector<int> SEQ_inrange(SEQ_combined.begin()+(reconstruction_start-StartPos),SEQ_combined.end());
								vector<int> Ns(gene_length-SEQ_inrange.size(),0);
								SEQ_inrange.insert(SEQ_inrange.end(),Ns.begin(),Ns.end());
								SEQ_range = SEQ_inrange;
							}
							else
							{
								vector<int> SEQ_inrange(SEQ_combined.begin()+(reconstruction_start-StartPos),SEQ_combined.begin()+(reconstruction_end-StartPos+1));
								SEQ_range = SEQ_inrange;
							}
						}
						else
						{
							if (EndPos <= reconstruction_end)
							{
								vector<int> SEQ_inrange = SEQ_combined;
								vector<int> Ns1(StartPos-reconstruction_start,0);
								SEQ_inrange.insert(SEQ_inrange.begin(),Ns1.begin(),Ns1.end());
								vector<int> Ns2(reconstruction_end-EndPos,0);
								SEQ_inrange.insert(SEQ_inrange.end(),Ns2.begin(),Ns2.end());	
								SEQ_range = SEQ_inrange;						
							}
							else
							{
								vector<int> SEQ_inrange(SEQ_combined.begin(),SEQ_combined.begin()+(reconstruction_end-StartPos+1));
								vector<int> Ns(StartPos-reconstruction_start,0);
								SEQ_inrange.insert(SEQ_inrange.begin(),Ns.begin(),Ns.end());
								SEQ_range = SEQ_inrange;
							}
						}
						Read_matrix.push_back(SEQ_range); // nReads by genome_length
			  			mean_length += SEQ_combined.size();
						seq_counter++;
						StartSeq.push_back(StartPos);//vector saving start positions for each paired-end sequence
                                                EndSeq.push_back(EndPos); //vector saving end positions for each paired-end sequence
                                                ReadSeq.push_back(SEQ_combined);//saving filtered paired-end sequences
					}				  
				}
	      		}
			else // unmapped or low qual/short seqlen etc
			{
				unmapped_counter++;
				if (qual < min_qual) // || seq_b.size() < min_length)
                                {
                                	vector<int> tag(8,0);
                                        tag[0] = pair_counter;
                                        tag[2] = qual;
					tag[3] = al_start;
					tag[4] = seq_b.size();
					tag[5] = qual_pairs;
					tag[6] = al_start_pairs;
					tag[7] = seq_b_pairs.size();
                                        lowQSseq.push_back(tag);
                                }
                                if (qual_pairs < min_qual) // || seq_b_pairs.size() < min_length)
                                {
                                        vector<int> tag(8,1);
                                        tag[0] = pair_counter;
                                        tag[2] = qual;
					tag[3] = al_start;
                                        tag[4] = seq_b.size();
                                        tag[5] = qual_pairs;
                                        tag[6] = al_start_pairs;
                                        tag[7] = seq_b_pairs.size();
                                        lowQSseq.push_back(tag);
                                }
			}
	   	 }
	}
	std::ofstream writefile1;
        std::ofstream writefile2;
        std::ofstream writefile3;
        std::ofstream writefile4;
	std::string name = zonename;
	writefile1.open(name+"_lowQSseq.txt");
	writefile2.open(name+"_StartSeq.txt");
        writefile3.open(name+"_EndSeq.txt");
        writefile4.open(name+"_ReadSeq.txt");
	for (int i = 0; i < lowQSseq.size(); i++)
	{
		writefile1 << lowQSseq[i][0] << " " << lowQSseq[i][1] << " " << lowQSseq[i][2] << " " << lowQSseq[i][3] << " " << lowQSseq[i][4] << " " << lowQSseq[i][5] << " " << lowQSseq[i][6] << " " << lowQSseq[i][7]  <<endl;
	}
	for (int i=0; i<seq_counter; i++)
		writefile2 << StartSeq[i] << " ";
	for (int i=0; i<seq_counter; i++)
                writefile3 << EndSeq[i] << " ";
	for (int i=0; i<seq_counter; i++)
	{
		for (int j=0; j<ReadSeq[i].size(); j++)
			writefile4 << ReadSeq[i][j] << " ";
		writefile4 << "\n";
	}
	writefile1.close();
	writefile2.close();
	writefile3.close();
	writefile4.close();
	cout << "mapped_counter:"<<mapped_counter<< " unmapped_counter:"<<unmapped_counter<< endl;
	cout << "num_lowQSseq:"<< lowQSseq.size()<<" seq_count:"<<seq_counter<<endl;
	cout << "filtered_counter1:"<< filtered_counter1 <<" filtered_counter2:"<<filtered_counter2<<endl;
        cout << " filtered_cond2:"<<filtered_cond2<<endl;
	cout << "pair_counter:" << pair_counter <<" singleton_counter:"<<singleton_counter<< " total_counter:"<<total_count<<endl;
  	return 0;
}

void callSNV(vector<vector<int> >& Read_matrix, vector<vector<int> >& Allele_freq, vector<vector<int> >& SNV_matrix,  vector<vector<int> >& SNV_freq, vector<int>& Homo_seq, vector<int>& SNV_pos, int nReads, int gene_length, double SNV_thres, int& nSNV, vector<int>& deleted_reads_list, string zonename)
{
	for (int j=0; j<gene_length; j++)
	{
		vector<int> count(5,0);
		for (int i=0; i<nReads; i++)
		{
			switch (Read_matrix[i][j])
			{
				case 1:
					count[0] = count[0] + 1;
					break;
				case 2:
					count[1] = count[1] + 1;
					break;
				case 3:
					count[2] = count[2] + 1;
					break;
	 			case 4:
					count[3] = count[3] + 1;
					break;			
			}
		}
		count[4] = count[0] + count[1] + count[2] +count[3];
		
		Allele_freq.push_back(count); // gene_length by 5
	}
	
	// get SNV position and corresponding allele_freq and SNV_matrix
	vector<vector<int> > temp_snv_matrix;
	for (int i = 0; i< gene_length; i++)
	{
		vector<int> allele_delete;
		int allele_sum = 0;
		int count = 0;
		if (Allele_freq[i][4] > 0)
		{
			for (int j = 0; j<4 ;j++)
			{
				double ratio = Allele_freq[i][j]/double(Allele_freq[i][4]);
				if (ratio < SNV_thres)
					allele_delete.push_back(0);
				else
				{
					allele_delete.push_back(Allele_freq[i][j]);
					allele_sum = allele_sum + Allele_freq[i][j];
					count++;
				}		
			}
			allele_delete.push_back(allele_sum);
		}
		// delete potential sequencing error
		vector<int> error_delete;
		if (count > 1)
		{
			SNV_freq.push_back(allele_delete); // nSNV by 5
			SNV_pos.push_back(i); // nSNV by 1
			
			for (int j = 0; j< nReads; j++)
			{
				if (Read_matrix[j][i]!=0)
				{
					if (allele_delete[Read_matrix[j][i]-1]!=0)
						error_delete.push_back(Read_matrix[j][i]);
					else
						error_delete.push_back(0);
				}
				else
					error_delete.push_back(0);
			}
			temp_snv_matrix.push_back(error_delete); // nSNV by nReads after error correction
			Homo_seq.push_back(0);
		}	
		else
		{
			int maxind = -1;
			int maxval = 0;
			for (int j = 0;j<4;j++)
			{
				if (maxval < Allele_freq[i][j])
				{
					maxval = Allele_freq[i][j];
					maxind = j;
				}
			}	
			Homo_seq.push_back(maxind+1);
		}
	}
	nSNV = SNV_pos.size();
	
	// delete empty fragment
	vector<int> ReadInd;
	for (int i=0; i< nReads; i++)
	{ 
		int count = 0;
		vector<int> frag_delete;
		for(int j=0; j<nSNV; j++)
		{
			if (temp_snv_matrix[j][i]!=0)
				count++;
		}
		if (count!=0)
		{	
			for(int j=0; j<nSNV; j++)
				frag_delete.push_back(temp_snv_matrix[j][i]);
			SNV_matrix.push_back(frag_delete);
			ReadInd.push_back(i+1);
		}
		else
			deleted_reads_list.push_back(i);
	}
	std::ofstream writefile1;
        std::string name = zonename;
        writefile1.open(name+"_ReadInd.txt");
	for (int i=0; i<ReadInd.size(); i++)
                writefile1 << ReadInd[i] << " ";
	writefile1.close();
}


int main(int argc, char* argv[]) {
	if(argc < 2)
	{
    		cout <<"usage: aBayesQR <config.txt>" << endl;
    		return 0;
  	}
 
  	srand(time(NULL));
 
  	string line, line_stats, line_ID;
  	string tok = ":";
  	vector<string> tokens, tokens2, tokens_as;
  
  	string confStr;
 
	if(argc == 2) confStr = argv[1];

	ifstream infConf(confStr.c_str(), ios::in);
  	vector<string> arg_buffer;
	int count =0;
	string arg_buffer_sub;
	
	cout << endl;
	while (!infConf.eof() )
	{
		getline(infConf,line,'\n');
	    	if( line.length() > 0 ) 
		{
	      		size_t pos = line.find(":");
	      		cout << count++ << ". "<< line << endl;
	      		arg_buffer_sub = line.substr(pos+2);
	      		arg_buffer.push_back(arg_buffer_sub);
	      		arg_buffer_sub.clear();
	    	}
  	}
  	infConf.close();
	
	string cons;
	vector<string> FASTAreads;
	double SNV_thres;
	double overlap_ratio_start, overlap_ratio_end;
	int  reconstruction_start,  reconstruction_end;
	double  min_qual;
	int min_length, max_insertln;
	double seq_err;
	double eta1;
	int initial_k;
	
	string zonename;
	
	if(count < 3){
	    cout <<"problem with config file...please check" << endl;
	    return 0;
	}
	cons =  arg_buffer[0];
	FASTAreads.push_back(arg_buffer[1]);
	SNV_thres = atof(arg_buffer[2].c_str());

        if(count > 3)
	{
		reconstruction_start = atoi(arg_buffer[3].c_str());
		if(count > 4)
		{
			reconstruction_end = atoi(arg_buffer[4].c_str());
		  	if(count > 5)
			{
		    		min_qual = atoi(arg_buffer[5].c_str());
		    		if(count > 6)
				{
		      			min_length =  atoi(arg_buffer[6].c_str());
		      			if(count > 7)
					{
		      				max_insertln = atoi(arg_buffer[7].c_str());
                        			if(count >8)
						{
                          				zonename = arg_buffer[8];
                          				if(count > 9)
							{
                            					seq_err = atof(arg_buffer[9].c_str());
                            					if(count > 10)
								{
									eta1 = atof(arg_buffer[10].c_str());
									if(count > 11)
									{
                              							initial_k = atof(arg_buffer[11].c_str());
									}
                            					}
                          				}
                        			}
		      			}
		    		}
		  	}
		}
	}
	seq_err = seq_err*0.01;	
	
	clock_t start, end;
	double cpu_time_used;
	start = clock();

	vector<vector<int> > Read_matrix;	  
	int total_count = 0;
	char gap_quality = '*'; //'I'; 	 
	double mean_length = 0.;
	int gene_length  = reconstruction_end-reconstruction_start+1;
        vector<vector<int> > ReadSeq;
        vector<int> StartSeq, EndSeq;	

	int error_flag = 0;
  	error_flag = parseSAMpaired(FASTAreads[0], min_qual,  min_length, max_insertln, gap_quality,  mean_length, Read_matrix, reconstruction_start, reconstruction_end, total_count, gene_length, zonename, ReadSeq, StartSeq, EndSeq);
 	if(error_flag>0)
    		return 1;
	int nReads = Read_matrix.size();

	for (int i = 0; i< Read_matrix.size(); i++)
	{
	        if (Read_matrix[i].size() != gene_length)
			cout << "error!!!"<<endl;
	}
	cout <<  endl << "After parsing " << total_count << " reads in file " <<FASTAreads[0]<< ", there are "<<nReads<< " paired_end reads(mean lengths "<< mean_length/Read_matrix.size() << ") covering regions "<< reconstruction_start << "-" << reconstruction_end <<"."<< endl;
	
	end = clock();	
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cout << "CPU time for SAM parsing: "  << cpu_time_used << endl<<endl;

	start = clock();		
	vector<vector<int> > Allele_freq; // nReads by 5
	vector<vector<int> > SNV_matrix; // nFrag by nSNV
	vector<vector<int> > SNV_freq; // nSNV by 5
	vector<int> Homo_seq;
	vector<int> SNV_pos;
	int nSNV;
	vector<int> deleted_reads_list;	
	callSNV(Read_matrix, Allele_freq, SNV_matrix, SNV_freq, Homo_seq, SNV_pos, nReads, gene_length, SNV_thres, nSNV,deleted_reads_list, zonename);
	int nFrag = SNV_matrix.size();	
	cout << "After calling SNVs from " << gene_length << " bases in regions between " << reconstruction_start << " and " << reconstruction_end << ", " << nSNV<< " SNVs are detected." << endl;
	cout << "After correcting error, "<< nFrag <<" fragments are used for quasi-species reconstruction." << endl;
	cout << "reduced number of fragment:" << deleted_reads_list.size() << endl;
	end = clock();	
	cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
	cout << "CPU time for SNV calling: "  << cpu_time_used << endl << endl;
	

	long long int read_ln = min_length;
        std::string name = zonename;

	std::ofstream writefile1;
	std::ofstream writefile2;
	std::ofstream writefile3;

	writefile1.open(name+"_SNV_pos.txt");
	writefile2.open(name+"_SNV_matrix.txt");
	writefile3.open(name+"_Homo_seq.txt");

	for(int i=0; i<SNV_pos.size();i++)
                writefile1<<SNV_pos[i]<<" ";
	cout << "save SNV_pos"<<endl;

	for(int i=0; i<nFrag;i++)
        {
                for (int j = 0; j<nSNV; j++)
                        writefile2 << SNV_matrix[i][j]<<" ";
                writefile2 << "\n";
        }
	cout << "save SNV_matrix"<<endl;
	
	for (int i=0; i<Homo_seq.size(); i++)
		writefile3<<Homo_seq[i]<<" ";
	cout << "save Home_seq"<< endl;
	
	writefile1.close();
	writefile2.close();
	writefile3.close();
}
