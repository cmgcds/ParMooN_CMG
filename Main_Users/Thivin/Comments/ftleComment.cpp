
   // // DOF MAPPING
    // cout << "DOF MAPPING " <<endl;
    // for( int i = 0 ; i < DOF2CellMapping.size() ; i++)
    // {
    //     cout << " DOF : " << i <<" --> ";
    //     for ( int j = 0 ; j < DOF2CellMapping[i].size(); j++)
    //     {
    //         cout <<  ", "<<DOF2CellMapping[i][j] ;
    //     }
    //     cout <<endl;
    // }

    // // CELL MAPPING
    // cout << "CELL MAPPING " <<endl;
    // for( int i = 0 ; i < Cell2DOFMapping.size() ; i++)
    // {
    //     cout << " CELL : " << i <<" --> " ;
    //     for ( int j = 0 ; j < Cell2DOFMapping[i].size(); j++)
    //     {
    //         cout <<  ", "<<Cell2DOFMapping[i][j] ;
    //     }
    //     cout <<endl;
    // }



// Code comment snippets for FTLE 

   for(int i = 0 ; i < N_Particles; i++)
    {
        cout << "--- Particle " << i << " ------- "<<endl;
        cout << " Pos : " << x_pos[i] << " , " << y_pos[i] <<endl;

        cout << " Left : " ;
        if(neibhourLeft[i] >= 0)
            cout << x_pos[neibhourLeft[i]] << " , " << y_pos[neibhourLeft[i]]  << " DOF : " << neibhourLeft[i] ;
        else
            cout << " *** Corner left *** ";
        cout <<endl;

        cout << " Right : " ;
        if(neibhourRight[i] >=0)
            cout << x_pos[neibhourRight[i]] << " , " << y_pos[neibhourRight[i]] << " DOF : " << neibhourRight[i];
        else
            cout << " *** Corner Right *** ";
        cout <<endl;

        cout << " Top : " ;
        if(neibhourTop[i] >=0)
            cout << x_pos[neibhourTop[i]] << " , " << y_pos[neibhourTop[i]] << " DOF : " << neibhourTop[i];
        else
            cout << " *** Corner Top *** ";
        cout <<endl;

        cout << " Bottom : " ;
        if(neibhourBottom[i]>=0)
            cout << x_pos[neibhourBottom[i]] << " , " << y_pos[neibhourBottom[i]]<< " DOF : " << neibhourBottom[i] ;
        else
            cout << " *** Corner Bottom *** ";
        cout <<endl;
    }