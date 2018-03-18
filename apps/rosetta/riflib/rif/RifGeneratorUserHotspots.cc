// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.



#include <riflib/rif/RifGeneratorUserHotspots.hh>


	#include <ObjexxFCL/format.hh>

	#include <boost/random/mersenne_twister.hpp>
	#include <boost/random/uniform_real.hpp>

	#include <core/id/AtomID.hh>
	#include <core/pose/Pose.hh>
    #include <core/import_pose/import_pose.hh>
	#include <core/scoring/motif/util.hh>

	#include <devel/init.hh>
	#include <riflib/RotamerGenerator.hh>
	#include <riflib/util.hh>

    #include <scheme/actor/Atom.hh>
    #include <scheme/numeric/rand_xform.hh>


namespace devel {
namespace scheme {
namespace rif {



	void
	RifGeneratorUserHotspots::generate_rif(
		RifAccumulatorP accumulator,
		RifGenParamsP params
	){

		typedef numeric::xyzVector<core::Real> Vec;

		// some sanity checks
		int const n_hspot_groups = this->opts.hotspot_files.size();
		runtime_assert_msg( n_hspot_groups, "no hotspot group files specified!!" );
		runtime_assert_msg( n_hspot_groups<16, "too many hotspot groups!!" );
		// runtime_assert_msg( accumulator->rif()->has_sat_data_slots(), "This RIF type doesn't support sat groups!!!" );

		std::cout << "RifGeneratorUserHotspots opts:" << std::endl;
		std::cout << "    hotspot_sample_cart_bound:  " << this->opts.hotspot_sample_cart_bound << std::endl;
        std::cout << "    hotspot_sample_angle_bound: " << this->opts.hotspot_sample_angle_bound << std::endl;
        std::cout << "    hotspot_nsamples:           " << this->opts.hotspot_nsamples << std::endl;
		std::cout << "    hbond_weight:               " << this->opts.hbond_weight << std::endl;
		std::cout << "    upweight_multi_hbond:       " << this->opts.upweight_multi_hbond << std::endl;
		std::cout << "    target_center:              "
			<< this->opts.target_center[0] << " "
			<< this->opts.target_center[1] << " "
			<< this->opts.target_center[2] << std::endl;
        numeric::xyzVector<double> xyz_tgt_cen( this->opts.target_center[0], this->opts.target_center[1], this->opts.target_center[2] );
		for( auto s : this->opts.hotspot_files ){
			std::cout << "    hotspot_group:              " << s << std::endl;
		}


		// setup the hacky but fast scorer
		devel::scheme::ScoreRotamerVsTarget<
				VoxelArrayPtr, ::scheme::chemical::HBondRay, ::devel::scheme::RotamerIndex
			> rot_tgt_scorer;
		{
			std::vector< ::scheme::chemical::HBondRay > target_donors, target_acceptors;
			for( auto ir : params->target_res ){
				::devel::scheme::get_donor_rays   ( *params->target, ir, params->hbopt, target_donors );
				::devel::scheme::get_acceptor_rays( *params->target, ir, params->hbopt, target_acceptors );
			}
			std::cout << "target_donors.size() " << target_donors.size() << " target_acceptors.size() " << target_acceptors.size() << std::endl;
			{
				rot_tgt_scorer.rot_index_p_ = params->rot_index_p;
				rot_tgt_scorer.target_field_by_atype_ = params->field_by_atype;
				rot_tgt_scorer.target_donors_ = target_donors;
				rot_tgt_scorer.target_acceptors_ = target_acceptors;
				rot_tgt_scorer.hbond_weight_ = this->opts.hbond_weight;
				rot_tgt_scorer.upweight_multi_hbond_ = this->opts.upweight_multi_hbond;
				rot_tgt_scorer.upweight_iface_ = 1.0;

			}
		}

		std::vector<EigenXform> sample_position_deltas{ EigenXform::Identity() };
        std::mt19937 rng((unsigned int)time(0) + 293754);
        for( int i = 0; i < this->opts.hotspot_nsamples-1; ++i ){
            EigenXform xrand;
            float ang_bound_radians = this->opts.hotspot_sample_angle_bound / 180.0 * M_PI;
            ::scheme::numeric::rand_xform_sphere(rng, xrand, this->opts.hotspot_sample_cart_bound, ang_bound_radians );
            sample_position_deltas.push_back(xrand);
        }
        // std::cout << sample_position_deltas.size() << std::endl;
        // utility_exit_with_message("foianrst");

        if( this->opts.dump_hotspot_samples ){
            params->target->dump_pdb("target.pdb");
        }


		// fill this in somehow... maybe start with random purterbations... later I can help you do it "right" with a grid of some kind

		// loop over files (one file is one hotspot group)
		for( int i_hotspot_group = 0; i_hotspot_group < this->opts.hotspot_files.size(); ++i_hotspot_group ){

            std::string const & hotspot_file = this->opts.hotspot_files[i_hotspot_group];
            core::pose::Pose pose;
            core::import_pose::pose_from_file( pose, hotspot_file );

            std::cout << "i_hotspot_group " << i_hotspot_group << " " << hotspot_file << std::endl;

			// read in pdb files # i_hotspot_group
			for( int i_hspot_res = 1; i_hspot_res <= pose.size(); ++i_hspot_res ){
                std::cout << "    i_hspot_res " << i_hspot_res << " " << hotspot_file << " res " << i_hspot_res << std::endl;

                int input_nheavy = pose.residue(i_hspot_res).nheavyatoms();
                std::cout << "    align atoms:"
                          << " " << pose.residue(i_hspot_res).atom_name(input_nheavy-2)
                          << " " << pose.residue(i_hspot_res).atom_name(input_nheavy-1)
                          << " " << pose.residue(i_hspot_res).atom_name(input_nheavy-0)
                          << std::endl;
                // "stub" to align to defined by last three heavy atoms
                EigenXform Xref = ::scheme::chemical::make_stub<EigenXform>(
                    pose.residue(i_hspot_res).xyz( input_nheavy - 2 ) - xyz_tgt_cen,
                    pose.residue(i_hspot_res).xyz( input_nheavy - 1 ) - xyz_tgt_cen,
                    pose.residue(i_hspot_res).xyz( input_nheavy - 0 ) - xyz_tgt_cen
                );

				// for each irot that is the right restype (can be had from rot_intex_p)
                int irot_begin = params->rot_index_p->index_bounds(pose.residue(i_hspot_res).name3()).first;
                int irot_end   = params->rot_index_p->index_bounds(pose.residue(i_hspot_res).name3()).second;
                if( irot_begin >= irot_end ){
                    std::cerr << "unknown residue " << pose.residue(i_hspot_res).name3() << std::endl;
                    std::exit(-1);
                }

				for( int irot = irot_begin; irot < irot_end; ++irot ){
                    std::cout << "        irot " << irot << std::endl;

					std::vector<SchemeAtom> const & rotamer_atoms( params->rot_index_p->atoms(irot) );
                    // "stub" to align to defined by last three heavy atoms, note 0-indexing
                    EigenXform Xrotamer = ::scheme::chemical::make_stub<EigenXform>(
                        rotamer_atoms.at( params->rot_index_p->nheavyatoms(irot) - 3 ).position(),
                        rotamer_atoms.at( params->rot_index_p->nheavyatoms(irot) - 2 ).position(),
                        rotamer_atoms.at( params->rot_index_p->nheavyatoms(irot) - 1 ).position()
                    );

					// figure out the transform that aligns the rotamer's standard position onto your input hotspot
					// res to align to will be pose.residue(i_hspot_res)
					// this is a dummy
					EigenXform x_orig_position = Xref * Xrotamer.inverse();

					for( int i_ptrb = 0; i_ptrb < sample_position_deltas.size(); ++i_ptrb ){

                        // TODO: fix to rotate around COM
						EigenXform x_position = x_orig_position * sample_position_deltas[i_ptrb];

						// you can check their "energies" against the target like this, obviously substituting the real rot# and position
						float positioned_rotamer_score = rot_tgt_scorer.score_rotamer_v_target( irot, x_position );


						// add the rotamer to the rif if it's any good
						if( positioned_rotamer_score <= this->opts.hotspot_score_thresh ){ // probably want this threshold to be an option or something
                            // std::cout << "            positioned_rotamer_score " << positioned_rotamer_score << std::endl;

                            if( this->opts.dump_hotspot_samples ){
                                std::string outfname = "test_hs_gen_"+utility::file_basename(hotspot_file)+"_"+str(i_hspot_res)+"_"+str(irot)+"_"+str(i_ptrb)+".pdb";
                                // std::cout << "dumping " << outfname << std::endl;
                                std::ofstream out(outfname);
                                out << "MODEL" << std::endl;
                                for( auto a : rotamer_atoms ){
                                    a.set_position( x_position * a.position() );
                                    ::scheme::actor::write_pdb( out, a, params->rot_index_p->chem_index_ );
                                }
                                out << "ENDMDL" << std::endl;
                                out.close();
                            }

							accumulator->insert( x_position, positioned_rotamer_score, irot, i_hotspot_group, -1 );


                            accumulator->checkpoint( std::cout );



						}

					} // end position perturbations

				} // end loop over rotamers which match hotspot res

			} //  end loop over residues in hotspot group

		} // end loop over hotspot groups

		// let the rif builder thing know you're done
        accumulator->checkpoint( std::cout );


        auto rif_ptr = accumulator->rif();
        std::cout << "testing rifbase key iteration" << std::endl;
        int count = 0;
        for( auto key : rif_ptr->key_range() ){
            EigenXform bin_center = rif_ptr->get_bin_center(key);
            // right... the edges of the bins.... this is only *mostly* true
            // runtime_assert( rif_ptr->get_bin_key(bin_center) == key );
            std::cout << "BIN: key " << key << " xform.trans: " << bin_center.translation().transpose() << std::endl;
            auto rotscores = rif_ptr->get_rotamers_for_key(key);
            runtime_assert( rotscores.size() > 0 );
            for( auto rot_score : rotscores ){
                // can't wait for cxx17 structured bindings!!!
                float rotamer_score = rot_score.first;
                int rotamer_number = rot_score.second;
                std::string resn = params->rot_index_p->resname(rotamer_number);
                std::cout << " rotamer " << rotamer_number << " " << resn << ", score " << rotamer_score << std::endl;
            }

            if(++count > 10) utility_exit_with_message("aireost");
        }


        // std::exit(-1);
	}


}
}
}

