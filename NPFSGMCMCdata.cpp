//
//  NPFSGMCMCdata.cpp
//  
//
//  Created by Ayan Acharya on 9/20/15.
//
//

#include "NPFSGMCMCdata.h"

data::data(string trFileName)
{
    string line;
    ifstream trFile(trFileName); 
            
    int count, d, w, p, q, tmp1, tmp2, tmp3, tmp4, newdoc;
    
    cout<<"loading of data starts for NPF-SGMCMC .."<<endl;
    
    getline(trFile, line); istringstream iss1(line); iss1 >> D >> V >> S; 
    Xdw = sp_mat(D,V); 
    
    //cout<<D<<"\t"<<V<<"\t"<<endl;    
    
    tmp1 = 0; tmp3 = -1; d = 0;
    while (getline(trFile, line))
    {
        std::istringstream iss1(line);
        if (!(iss1 >> tmp1 >> tmp2 >> tmp4))
        {
            break;
        }
        if(tmp1>tmp3)
        {
            tmp3 = tmp1;
            if(tmp1>0)
				d += 1;
			Xdw(d,tmp2) = tmp4;     
        }
        else
            Xdw(d,tmp2) = tmp4;
    }
    
    cout<<"loading of data ends for NPF-SGMCMC .."<<endl;
    return;
};
