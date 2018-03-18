#include <basic/options/option_macros.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <vector>


OPT_1GRP_KEY(     StringVector , rif_dock, scaffolds )
	OPT_1GRP_KEY(  StringVector, rif_dock, scaffold_res )
	OPT_1GRP_KEY(  StringVector, rif_dock, scaffold_res_fixed )
	OPT_1GRP_KEY(  Boolean     , rif_dock, scaffold_res_use_best_guess )
	OPT_1GRP_KEY(  Boolean     , rif_dock, scaffold_to_ala )
	OPT_1GRP_KEY(  Boolean     , rif_dock, scaffold_to_ala_selonly )
	OPT_1GRP_KEY(  Boolean     , rif_dock, replace_orig_scaffold_res )
	OPT_1GRP_KEY(  Boolean     , rif_dock, replace_all_with_ala_1bre )
	OPT_1GRP_KEY(  Boolean     , rif_dock, random_perturb_scaffold )

	OPT_1GRP_KEY(  StringVector, rif_dock, target_bounding_xmaps )
	OPT_1GRP_KEY(  String      , rif_dock, target_pdb )
	OPT_1GRP_KEY(  String      , rif_dock, target_res )
	OPT_1GRP_KEY(  String      , rif_dock, target_rif )
	OPT_1GRP_KEY(  Real        , rif_dock, target_rf_resl )
	OPT_1GRP_KEY(  Integer     , rif_dock, target_rf_oversample )
	OPT_1GRP_KEY(  String      , rif_dock, target_rf_cache )

	OPT_1GRP_KEY(  StringVector, rif_dock, data_cache_dir )

	OPT_1GRP_KEY(  Real        , rif_dock, beam_size_M )
	OPT_1GRP_KEY(  Real        , rif_dock, search_diameter )
	OPT_1GRP_KEY(  Real        , rif_dock, hsearch_scale_factor )

	OPT_1GRP_KEY(  Real        , rif_dock, max_rf_bounding_ratio )
	OPT_1GRP_KEY(  Boolean     , rif_dock, make_bounding_plot_data )
	OPT_1GRP_KEY(  Boolean     , rif_dock, align_output_to_scaffold )
	OPT_1GRP_KEY(  Boolean     , rif_dock, output_scaffold_only )
	OPT_1GRP_KEY(  Integer     , rif_dock, n_pdb_out )

	OPT_1GRP_KEY(  Real        , rif_dock, rf_resl )
	OPT_1GRP_KEY(  Integer     , rif_dock, rf_oversample )
	OPT_1GRP_KEY(  Boolean     , rif_dock, downscale_atr_by_hierarchy )

	OPT_1GRP_KEY(  Integer     , rif_dock, rotrf_oversample )
	OPT_1GRP_KEY(  Real        , rif_dock, rotrf_resl )
	OPT_1GRP_KEY(  Real        , rif_dock, rotrf_spread )
	OPT_1GRP_KEY(  Real        , rif_dock, rotrf_scale_atr )
	OPT_1GRP_KEY(  String      , rif_dock, rotrf_cache_dir )

	OPT_1GRP_KEY(  Boolean     , rif_dock, hack_pack )
	OPT_1GRP_KEY(  Real        , rif_dock, hack_pack_frac )
	OPT_1GRP_KEY(  Real        , rif_dock, pack_iter_mult )
	OPT_1GRP_KEY(  Integer     , rif_dock, pack_n_iters )
	OPT_1GRP_KEY(  Real        , rif_dock, hbond_weight )
	OPT_1GRP_KEY(  Real        , rif_dock, upweight_multi_hbond )
	OPT_1GRP_KEY(  Real        , rif_dock, global_score_cut )

	OPT_1GRP_KEY(  Integer     , rif_dock, n_result_limit )
	OPT_1GRP_KEY(  Real        , rif_dock, redundancy_filter_mag )

	OPT_1GRP_KEY(  Real        , rif_dock, force_output_if_close_to_input )
	OPT_1GRP_KEY(  Integer     , rif_dock, force_output_if_close_to_input_num )

	OPT_1GRP_KEY(  Real        , rif_dock, upweight_iface )

	OPT_1GRP_KEY(  Boolean     , rif_dock, use_scaffold_bounding_grids )

	OPT_1GRP_KEY(  Boolean     , rif_dock, restrict_to_native_scaffold_res )
	OPT_1GRP_KEY(  Real        , rif_dock, bonus_to_native_scaffold_res )
	OPT_1GRP_KEY(  Boolean     , rif_dock, add_native_scaffold_rots_when_packing )

	OPT_1GRP_KEY(  Boolean     , rif_dock, dump_all_rif_rots )

	OPT_1GRP_KEY(  String     , rif_dock, dokfile )
	OPT_1GRP_KEY(  String     , rif_dock, outdir )
	OPT_1GRP_KEY(  String     , rif_dock, output_tag )

	OPT_1GRP_KEY(  Boolean    , rif_dock, dont_use_scaffold_loops )

	OPT_1GRP_KEY(  Boolean    , rif_dock, full_scaffold_output )
	OPT_1GRP_KEY(  Boolean    , rif_dock, dump_resfile )
	OPT_1GRP_KEY(  Boolean    , rif_dock, pdb_info_pikaa )

	OPT_1GRP_KEY(  Boolean    , rif_dock, cache_scaffold_data )

	OPT_1GRP_KEY(  Real        , rif_dock, tether_to_input_position )

	OPT_1GRP_KEY(  Boolean     , rif_dock, lowres_sterics_cbonly )

	OPT_1GRP_KEY(  Integer     , rif_dock, require_satisfaction )

	OPT_1GRP_KEY(  Real        , rif_dock, rosetta_score_fraction )
	OPT_1GRP_KEY(  Real        , rif_dock, rosetta_score_then_min_below_thresh )
	OPT_1GRP_KEY(  Integer     , rif_dock, rosetta_score_at_least )
	OPT_1GRP_KEY(  Integer     , rif_dock, rosetta_score_at_most )
	OPT_1GRP_KEY(  Real        , rif_dock, rosetta_min_fraction )
	OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_min_fix_target )
	OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_min_targetbb )
	OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_min_scaffoldbb )
	OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_min_allbb )
	OPT_1GRP_KEY(  Real        , rif_dock, rosetta_score_cut )
	OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_hard_min )
	OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_score_total )
	OPT_1GRP_KEY(  Boolean     , rif_dock, rosetta_score_ddg_only )
	OPT_1GRP_KEY(  Real        , rif_dock, rosetta_score_rifres_rifres_weight )
	OPT_1GRP_KEY(  Real        , rif_dock, rosetta_score_rifres_scaffold_weight )
	OPT_1GRP_KEY(  String      , rif_dock, rosetta_soft_score )
	OPT_1GRP_KEY(  String      , rif_dock, rosetta_hard_score )

	OPT_1GRP_KEY(  Boolean     , rif_dock, extra_rotamers )
	OPT_1GRP_KEY(  Boolean     , rif_dock, extra_rif_rotamers )
	OPT_1GRP_KEY(  Integer     , rif_dock, always_available_rotamers_level )
    OPT_1GRP_KEY(  Boolean     , rif_dock, packing_use_rif_rotamers )

    OPT_1GRP_KEY(  Integer     , rif_dock, nfold_symmetry )
    OPT_1GRP_KEY(  RealVector  , rif_dock, symmetry_axis )

    // constrain file 
	OPT_1GRP_KEY(  StringVector, rif_dock, cst_files )
	OPT_1GRP_KEY(  StringVector, rif_dock, seeding_files )

	void register_options() {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		NEW_OPT(  rif_dock::scaffolds, "" , utility::vector1<std::string>() );
		NEW_OPT(  rif_dock::scaffold_res, "" , utility::vector1<std::string>() );
		NEW_OPT(  rif_dock::scaffold_res_fixed, "" , utility::vector1<std::string>() );
		NEW_OPT(  rif_dock::scaffold_res_use_best_guess, "" , false );
		NEW_OPT(  rif_dock::scaffold_to_ala, "" , false );
		NEW_OPT(  rif_dock::scaffold_to_ala_selonly, "" , true );
		NEW_OPT(  rif_dock::replace_orig_scaffold_res, "", true );
		NEW_OPT(  rif_dock::replace_all_with_ala_1bre, "" , false );
		NEW_OPT(  rif_dock::random_perturb_scaffold, "" , false );

		NEW_OPT(  rif_dock::target_bounding_xmaps, "" , utility::vector1<std::string>() );
		NEW_OPT(  rif_dock::target_pdb, "" , "" );
		NEW_OPT(  rif_dock::target_res, "" , "" );
		NEW_OPT(  rif_dock::target_rif, "" , "" );
		NEW_OPT(  rif_dock::target_rf_resl, ""       , 0.25 );
		NEW_OPT(  rif_dock::target_rf_oversample, "" , 2 );
		NEW_OPT(  rif_dock::downscale_atr_by_hierarchy, "" , true );

		NEW_OPT(  rif_dock::target_rf_cache, "" , "NO_CACHE_SPECIFIED_ON_COMMAND_LINE" );

		NEW_OPT(  rif_dock::data_cache_dir, "" , utility::vector1<std::string>(1,"./") );
		NEW_OPT(  rif_dock::beam_size_M, "" , 10.000000 );
		NEW_OPT(  rif_dock::max_rf_bounding_ratio, "" , 4 );
		NEW_OPT(  rif_dock::make_bounding_plot_data, "" , false );
		NEW_OPT(  rif_dock::align_output_to_scaffold, "" , false );
		NEW_OPT(  rif_dock::output_scaffold_only, "" , false );
		NEW_OPT(  rif_dock::n_pdb_out, "" , 10 );

		NEW_OPT(  rif_dock::rf_resl, ""       , 0.25 );
		NEW_OPT(  rif_dock::rf_oversample, "" , 2 );

		NEW_OPT(  rif_dock::rotrf_oversample, "" , 2 );
		NEW_OPT(  rif_dock::rotrf_resl, "" , 0.3 );
		NEW_OPT(  rif_dock::rotrf_spread, "" , 0.0 );
		NEW_OPT(  rif_dock::rotrf_scale_atr, "" , 1.0 );
		NEW_OPT(  rif_dock::rotrf_cache_dir, "" , "./" );

		NEW_OPT(  rif_dock::hack_pack, "" , true );
		NEW_OPT(  rif_dock::hack_pack_frac, "" , 0.2 );
		NEW_OPT(  rif_dock::pack_iter_mult, "" , 2.0 );
		NEW_OPT(  rif_dock::pack_n_iters, "" , 1 );
		NEW_OPT(  rif_dock::hbond_weight, "" , 2.0 );
		NEW_OPT(  rif_dock::upweight_multi_hbond, "" , 0.0 );
		NEW_OPT(  rif_dock::global_score_cut, "" , 0.0 );

		NEW_OPT(  rif_dock::n_result_limit, "" , 2000000000 );

		NEW_OPT(  rif_dock::redundancy_filter_mag, "" , 1.0 );

		NEW_OPT(  rif_dock::force_output_if_close_to_input, "" , 1.0 );
		NEW_OPT(  rif_dock::force_output_if_close_to_input_num, "" , 0 );

		NEW_OPT(  rif_dock::upweight_iface, "", 1.2 );

		NEW_OPT(  rif_dock::use_scaffold_bounding_grids, "", false );

		NEW_OPT(  rif_dock::search_diameter, "", 150.0 );
		NEW_OPT(  rif_dock::hsearch_scale_factor, "global scaling of rotation/translation search grid", 1.0 );

		NEW_OPT(  rif_dock::restrict_to_native_scaffold_res, "aka structure prediction CHEAT", false );
		NEW_OPT(  rif_dock::bonus_to_native_scaffold_res, "aka favor native CHEAT", -0.3 );
		NEW_OPT(  rif_dock::add_native_scaffold_rots_when_packing, "CHEAT", false );

		NEW_OPT(  rif_dock::dump_all_rif_rots, "", false );

		NEW_OPT(  rif_dock::dokfile, "", "default.dok" );
		NEW_OPT(  rif_dock::outdir, "", "./" );
		NEW_OPT(  rif_dock::output_tag, "", "" );

		NEW_OPT(  rif_dock::dont_use_scaffold_loops, "", false );

		NEW_OPT(  rif_dock::full_scaffold_output, "", false );
		NEW_OPT(  rif_dock::dump_resfile, "", false );
		NEW_OPT(  rif_dock::pdb_info_pikaa, "", false );

		NEW_OPT(  rif_dock::cache_scaffold_data, "", false );

		NEW_OPT(  rif_dock::tether_to_input_position, "", -1.0 );

		NEW_OPT(  rif_dock::lowres_sterics_cbonly, "", true );

		NEW_OPT(  rif_dock::require_satisfaction, "", 0 );

		NEW_OPT(  rif_dock::rosetta_score_fraction  , "",  0.00 );
		NEW_OPT(  rif_dock::rosetta_score_then_min_below_thresh, "", -9e9 );
		NEW_OPT(  rif_dock::rosetta_score_at_least, "", -1 );
		NEW_OPT(  rif_dock::rosetta_score_at_most, "", 999999999 );
		NEW_OPT(  rif_dock::rosetta_min_fraction  , "",  0.1 );
		NEW_OPT(  rif_dock::rosetta_min_targetbb  , "",  false );
		NEW_OPT(  rif_dock::rosetta_min_scaffoldbb  , "",  false );
		NEW_OPT(  rif_dock::rosetta_min_allbb  , "",  false );
		NEW_OPT(  rif_dock::rosetta_min_fix_target, "",  false );
		NEW_OPT(  rif_dock::rosetta_score_cut  , "", -10.0 );
		NEW_OPT(  rif_dock::rosetta_hard_min  , "", false );
		NEW_OPT(  rif_dock::rosetta_score_total  , "", false );
		NEW_OPT(  rif_dock::rosetta_score_ddg_only  , "", false );
		NEW_OPT(  rif_dock::rosetta_score_rifres_rifres_weight, "", 0.75 );
		NEW_OPT(  rif_dock::rosetta_score_rifres_scaffold_weight, "", 0.5 );
		NEW_OPT(  rif_dock::rosetta_soft_score, "", "beta_soft" );
		NEW_OPT(  rif_dock::rosetta_hard_score, "", "beta" );

		NEW_OPT(  rif_dock::extra_rotamers, "", true );
		NEW_OPT(  rif_dock::extra_rif_rotamers, "", true );
		NEW_OPT(  rif_dock::always_available_rotamers_level, "", 0 );
        NEW_OPT(  rif_dock::packing_use_rif_rotamers, "", true );

        NEW_OPT(  rif_dock::nfold_symmetry, "", 1 );
        NEW_OPT(  rif_dock::symmetry_axis, "", utility::vector1<double>() );

        // constrain file names
		NEW_OPT(  rif_dock::cst_files, "" , utility::vector1<std::string>() );
		NEW_OPT(  rif_dock::seeding_files, "" , utility::vector1<std::string>() );

	}

struct RifDockOpt
{
	std::vector<std::string> scaffold_fnames;
	std::vector<std::string> scaffold_res_fnames;
	std::vector<std::string> data_cache_path;
	std::vector<std::string> rif_files;

    // constrain fiile names
	std::vector<std::string> cst_fnames;
    std::vector<std::string> seeding_fnames;

	bool        VERBOSE                              ;
	double      resl0                                ;
	int64_t     DIM                                  ;
	int64_t     DIMPOW2                              ;
	int64_t     beam_size                            ;
	bool        replace_all_with_ala_1bre            ;
	bool        lowres_sterics_cbonly                ;
	float       tether_to_input_position_cut         ;
	bool        tether_to_input_position             ;
	float       global_score_cut                     ;
	std::string target_pdb                           ;
	std::string outdir                               ;
	std::string output_tag                           ;
	std::string dokfile_fname                        ;
	bool        dump_all_rif_rots                    ;
	bool        add_native_scaffold_rots_when_packing;
	bool        restrict_to_native_scaffold_res      ;
	float       bonus_to_native_scaffold_res         ;
	float       hack_pack_frac                       ;
	float       hsearch_scale_factor                 ;
	float       search_diameter                      ;
	bool        use_scaffold_bounding_grids          ;
	bool        scaffold_res_use_best_guess          ;
	bool        scaff2ala                            ;
	bool        scaff2alaselonly                     ;
	bool        replace_orig_scaffold_res            ;
	int         require_satisfaction                 ;
	float       target_rf_resl                       ;
	bool        align_to_scaffold                    ;
	bool        output_scaffold_only                 ;
	bool        pdb_info_pikaa                       ;
	bool        full_scaffold_output                 ;
	bool        dump_resfile                         ;
	std::string target_res_fname                     ;
	int         target_rf_oversample                 ;
	float       max_rf_bounding_ratio                ;
	std::string target_rf_cache                      ;
	bool        downscale_atr_by_hierarchy           ;
	bool        random_perturb_scaffold              ;
	bool        dont_use_scaffold_loops              ;
	bool        cache_scaffold_data                  ;
	float       rf_resl                              ;
	bool        hack_pack                            ;
	int         rf_oversample                        ;

	int         rotrf_oversample                     ;
	float       rotrf_resl                           ;
	float       rotrf_spread                         ;
	std::string rotrf_cache_dir                      ;
	float       rotrf_scale_atr                      ;

	float       pack_iter_mult                       ;
	int         pack_n_iters                         ;
	float       hbond_weight                         ;
	float       upweight_iface                       ;
	float       upweight_multi_hbond                 ;
	int         n_result_limit                       ;
	float       redundancy_filter_mag                ;
	int         force_output_if_close_to_input_num   ;
	float       force_output_if_close_to_input       ;
	int         n_pdb_out                            ;
	bool        extra_rotamers                       ;
	bool        extra_rif_rotamers                   ;
	int         always_available_rotamers_level      ;
	int         packing_use_rif_rotamers             ;

	float       rosetta_score_fraction               ;
	float       rosetta_score_then_min_below_thresh  ;
	float       rosetta_score_at_least               ;
	float       rosetta_score_at_most                ;
	float       rosetta_min_fraction                 ;
	bool        rosetta_min_fix_target               ;
	bool        rosetta_min_targetbb                 ;
	bool        rosetta_min_scaffoldbb               ;
	bool        rosetta_min_allbb                    ;
	float       rosetta_score_cut                    ;
	float       rosetta_hard_min                     ;
	bool        rosetta_score_total                  ;
	bool        rosetta_score_ddg_only               ;
	float       rosetta_score_rifres_rifres_weight   ;
	float       rosetta_score_rifres_scaffold_weight ;

	bool        rosetta_beta                         ;
	std::string rosetta_soft_score                   ;
	std::string rosetta_hard_score                   ;

    int         nfold_symmetry                       ;
    std::vector<float> symmetry_axis                 ;


	void init_from_cli()
	{
		using basic::options::option;
		using namespace basic::options::OptionKeys;

		runtime_assert( option[rif_dock::target_rif].user() );

		VERBOSE                                = false;
		resl0                                  = 16.0;
		DIM                                    = 6;
		DIMPOW2                                = 1<<DIM;
		beam_size                              = int64_t( option[rif_dock::beam_size_M]() * 1000000.0 / DIMPOW2 ) * DIMPOW2;
		replace_all_with_ala_1bre              = option[rif_dock::replace_all_with_ala_1bre          ]();

		target_pdb                             = option[rif_dock::target_pdb                         ]();
		lowres_sterics_cbonly                  = option[rif_dock::lowres_sterics_cbonly              ]();
		tether_to_input_position_cut           = option[rif_dock::tether_to_input_position           ]();
		tether_to_input_position               = tether_to_input_position_cut > 0.0;
		global_score_cut                       = option[rif_dock::global_score_cut                   ]();
		outdir                                 = option[rif_dock::outdir                             ]();
		output_tag                             = option[rif_dock::output_tag                         ]();
		dokfile_fname                          = outdir + "/" + option[rif_dock::dokfile             ]();
		dump_all_rif_rots                      = option[rif_dock::dump_all_rif_rots                  ]();
		add_native_scaffold_rots_when_packing  = option[rif_dock::add_native_scaffold_rots_when_packing ]();
		restrict_to_native_scaffold_res        = option[rif_dock::restrict_to_native_scaffold_res       ]();
		bonus_to_native_scaffold_res           = option[rif_dock::bonus_to_native_scaffold_res          ]();
		hack_pack_frac                         = option[rif_dock::hack_pack_frac                        ]();
		hsearch_scale_factor                   = option[rif_dock::hsearch_scale_factor                  ]();
		search_diameter                        = option[rif_dock::search_diameter                       ]();
		use_scaffold_bounding_grids            = option[rif_dock::use_scaffold_bounding_grids           ]();
		scaffold_res_use_best_guess            = option[rif_dock::scaffold_res_use_best_guess           ]();
		scaff2ala                              = option[rif_dock::scaffold_to_ala                       ]();
		scaff2alaselonly                       = option[rif_dock::scaffold_to_ala_selonly               ]();
		replace_orig_scaffold_res              = option[rif_dock::replace_orig_scaffold_res             ]();
		require_satisfaction                   = option[rif_dock::require_satisfaction                  ]();
		target_rf_resl                         = option[rif_dock::target_rf_resl                        ]();
		align_to_scaffold                      = option[rif_dock::align_output_to_scaffold              ]();
		output_scaffold_only                   = option[rif_dock::output_scaffold_only                  ]();
		pdb_info_pikaa                         = option[rif_dock::pdb_info_pikaa                        ]();
		full_scaffold_output                   = option[rif_dock::full_scaffold_output                  ]();
		dump_resfile                           = option[rif_dock::dump_resfile                          ]();
		target_res_fname                       = option[rif_dock::target_res                            ]();
		target_rf_oversample                   = option[rif_dock::target_rf_oversample                  ]();
		max_rf_bounding_ratio                  = option[rif_dock::max_rf_bounding_ratio                 ]();
		target_rf_cache                        = option[rif_dock::target_rf_cache                       ]();
		downscale_atr_by_hierarchy             = option[rif_dock::downscale_atr_by_hierarchy            ]();
		random_perturb_scaffold                = option[rif_dock::random_perturb_scaffold               ]();
		dont_use_scaffold_loops                = option[rif_dock::dont_use_scaffold_loops               ]();
		cache_scaffold_data                    = option[rif_dock::cache_scaffold_data                   ]();
		rf_resl                                = option[rif_dock::rf_resl                               ]();
		hack_pack                              = option[rif_dock::hack_pack                             ]();
		rf_oversample                          = option[rif_dock::rf_oversample                         ]();
		redundancy_filter_mag                  = option[rif_dock::redundancy_filter_mag                 ]();
		rotrf_oversample                       = option[rif_dock::rotrf_oversample                      ]();
		rotrf_resl                             = option[rif_dock::rotrf_resl                            ]();
		rotrf_spread                           = option[rif_dock::rotrf_spread                          ]();
		rotrf_cache_dir                        = option[rif_dock::rotrf_cache_dir                       ]();
		rotrf_scale_atr                        = option[rif_dock::rotrf_scale_atr                       ]();
		pack_iter_mult                         = option[rif_dock::pack_iter_mult                        ]();
		pack_n_iters                           = option[rif_dock::pack_n_iters                         ]();
		hbond_weight                           = option[rif_dock::hbond_weight                          ]();
		upweight_iface                         = option[rif_dock::upweight_iface                        ]();
		upweight_multi_hbond                   = option[rif_dock::upweight_multi_hbond                  ]();
		n_result_limit                         = option[rif_dock::n_result_limit                        ]();
		redundancy_filter_mag                  = option[rif_dock::redundancy_filter_mag                 ]();
		force_output_if_close_to_input_num     = option[rif_dock::force_output_if_close_to_input_num    ]();
		force_output_if_close_to_input         = option[rif_dock::force_output_if_close_to_input        ]();
		n_pdb_out                              = option[rif_dock::n_pdb_out                             ]();
		extra_rotamers                         = option[rif_dock::extra_rotamers                        ]();
		extra_rif_rotamers                     = option[rif_dock::extra_rif_rotamers                    ]();
		always_available_rotamers_level        = option[rif_dock::always_available_rotamers_level       ]();
		packing_use_rif_rotamers               = option[rif_dock::packing_use_rif_rotamers              ]();

  		rosetta_score_fraction                 = option[rif_dock::rosetta_score_fraction                ]();
  		rosetta_score_then_min_below_thresh    = option[rif_dock::rosetta_score_then_min_below_thresh   ]();
  		rosetta_score_at_least                 = option[rif_dock::rosetta_score_at_least                ]();
  		rosetta_score_at_most                  = option[rif_dock::rosetta_score_at_most                 ]();
  		rosetta_min_fraction                   = option[rif_dock::rosetta_min_fraction                  ]();
  		rosetta_min_fix_target                 = option[rif_dock::rosetta_min_fix_target                ]();
  		rosetta_min_targetbb                   = option[rif_dock::rosetta_min_targetbb                  ]();
  		rosetta_min_scaffoldbb                 = option[rif_dock::rosetta_min_scaffoldbb                ]();
  		rosetta_min_allbb                      = option[rif_dock::rosetta_min_allbb                     ]();
  		rosetta_score_cut                      = option[rif_dock::rosetta_score_cut                     ]();
  		rosetta_hard_min                       = option[rif_dock::rosetta_hard_min                      ]();
  		rosetta_score_total                    = option[rif_dock::rosetta_score_total                   ]();
  		rosetta_score_ddg_only                 = option[rif_dock::rosetta_score_ddg_only                ]();
  		rosetta_score_rifres_rifres_weight     = option[rif_dock::rosetta_score_rifres_rifres_weight    ]();
		rosetta_score_rifres_scaffold_weight   = option[rif_dock::rosetta_score_rifres_scaffold_weight  ]();
		rosetta_soft_score                     = option[rif_dock::rosetta_soft_score  ]();
		rosetta_hard_score                     = option[rif_dock::rosetta_hard_score  ]();
		rosetta_beta                           = option[corrections::beta]();


		for( std::string s : option[rif_dock::scaffolds     ]() )     scaffold_fnames.push_back(s);
		for( std::string s : option[rif_dock::scaffold_res  ]() ) scaffold_res_fnames.push_back(s);
		for( std::string s : option[rif_dock::data_cache_dir]() )     data_cache_path.push_back(s);

        // constrain file names
		for( std::string s : option[rif_dock::cst_files  ]() ) cst_fnames.push_back(s);
		for( std::string s : option[rif_dock::seeding_files  ]() ) seeding_fnames.push_back(s);

		for( std::string fn : option[rif_dock::target_bounding_xmaps]() ) rif_files.push_back(fn);
		rif_files.push_back( option[rif_dock::target_rif]() );

		if( scaff2ala && scaff2alaselonly &&  option[rif_dock::scaffold_to_ala_selonly].user() ){
			std::cout << "WARNING: -scaffold_to_ala overrides -scaffold_to_ala_selonly!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		}

		if( rosetta_score_total && rosetta_score_ddg_only ){
			std::cout << "WARNING: rosetta_score_total overrives rosetta_score_ddg_only" << std::endl;
			rosetta_score_ddg_only = false;
		}

        nfold_symmetry = option[rif_dock::nfold_symmetry]();
        symmetry_axis.clear();
        if( option[rif_dock::symmetry_axis]().size() == 3 ){
            symmetry_axis.push_back( option[rif_dock::symmetry_axis]()[1] );
            symmetry_axis.push_back( option[rif_dock::symmetry_axis]()[2] );
            symmetry_axis.push_back( option[rif_dock::symmetry_axis]()[3] );
        } else if( option[rif_dock::symmetry_axis]().size() == 0 ){
            symmetry_axis.push_back(0);
            symmetry_axis.push_back(0);
            symmetry_axis.push_back(1);
        } else {
            std::cout << "bad rif_dock::symmetry_axis option" << std::endl;
            std::exit(-1);
        }

	}
};

