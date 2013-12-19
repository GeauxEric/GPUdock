	for (conf = begin_row; conf <= end_row; ++conf) {

		DEBUG_4_("rank: ", rank);
		DEBUG_4_("conf", conf);
		n1 = ((atoms + 3) * conf);

		newpens = atoi((data[n1]).c_str());
		newlens = atoi((data[n1 + 1]).c_str());

#if TUNE
		d_complex->setProteinEnsembleCurrent(newpens);
		d_complex->setLigandEnsembleCurrent(newlens);
#endif

		for (int j = n1 + 2, k = 0; k < atoms; j++, k++) {
			stringstream vec1_entry((data[j]).c_str());
			vec1_entry >> newcoords[k][0] >> newcoords[k][1] >> newcoords[k][2];
		}

#if TUNE
		d_complex->setConfCoords(newcoords);	// set lig conf coords
		d_complex->clearContactMatrix();
		d_complex->createContactMatrix();
#endif

		int pens_p;
		int lens_p;
		Complex complex_p = *d_complex;
		double coords_p[MAXLIG][3];

		int conf2_p;
		int n;		//used for vector entries

		for (conf2_p = conf; conf2_p < total_confs; conf2_p++) {

#if 0
			DEBUG_3_("rank:", rank);
			DEBUG_3_("cycle: ", cycle);
			DEBUG_3_("conf2_p: ", conf2_p);
#endif

			if (conf2_p < total_confs) {

				//Calculate mcc for each conformation relative the the current "native"

				n = ((atoms + 3) * conf2_p);

				pens_p = atoi((data[n]).c_str());
				lens_p = atoi((data[n + 1]).c_str());

#if TUNE
				complex_p.setProteinEnsembleCurrent(pens_p);
				complex_p.setLigandEnsembleCurrent(lens_p);
#endif 

				for (int j = n + 2, k = 0; k < atoms; j++, k++) {
					stringstream vec_entry((data[j]).c_str());
					vec_entry >> coords_p[k][0] >> coords_p[k][1] >> coords_p[k][2];
				}

#if TUNE
				complex_p.setConfCoords(coords_p);
				complex_p.calculateEnergy();
#endif

				int mcc_index = 0;

				if (conf == 0) {
					mcc_index = conf2_p;	// does rank start at 0 or 1?
				} else {
					for (int i = 1; i <= conf; i++) {
						mcc_index += total_confs - (i - 1);
					}
					mcc_index += conf2_p - conf;
				}

				// cout <<  rank << conf << conf2_p << " mcc_index: " << mcc_index <<" cmcc " << complex_p.getEnergy(1) << endl;

#if 0
				DEBUG_4_("conf: ", conf);
				DEBUG_4_("mcc_index: ", mcc_index);
				DEBUG_4_("cmcc: ", complex_p.getEnergy(1));
				DEBUG_4_("rank: ", rank);
				DEBUG_4_("conf2_p: ", conf2_p);
#endif

				mccMatrix[mcc_index] = complex_p.getEnergy(1);
			}
		}


		// compute_a_row (row, mcc);
		// compute_a_row (total_confs - row - 1, mcc);
	}
