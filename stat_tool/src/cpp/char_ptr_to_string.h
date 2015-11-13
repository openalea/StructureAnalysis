#include <stat_tool/compound.h>
#include <stat_tool/multivariate_mixture.h>
#include <stat_tool/regression.h>
#include <stat_tool/convolution.h>
#include <stat_tool/stat_tools.h>
#include <stat_tool/curves.h>
#include <stat_tool/mixture.h>

namespace stat_tool
{
	bool spreadsheet_write_5245e9312629559f86ec2be0f532879d(class ::stat_tool::StatInterface const & self_57e8416c0462569e8916250db94d4d5a, class ::stat_tool::StatError & error, std::string const & path);
	
	bool ascii_write_235b9507a5a350ecab65344434d5cc94(class ::stat_tool::MultivariateMixtureData const & self_5cde05a9631d5dcebc04732a82e9a74e, class ::stat_tool::StatError & error, std::string const & path, bool exhaustive);
	
	bool plot_write_7b9e71cba8e25ee78598b85db88df807(class ::stat_tool::DiscreteDistributionData const & self_699cdb7e19d058c5bb8d30066fdc55ac, class ::stat_tool::StatError & error, std::string const & prefix, std::string const & title);
	
	void correction_update_b8755eed19cf514f8f6fe05624e5c415(class ::stat_tool::StatError & self_4cdf3268c85b531ab98dd99593a41730, std::string const & ilabel, std::string const & correction, int iline, int icolumn);
	
	bool plot_write_b3d8fdf8d47f52c5b8ddd7339945d686(class ::stat_tool::Mixture const & self_7fd12c45b7445cbda440f447ba16a244, class ::stat_tool::StatError & error, std::string const & prefix, std::string const & title);
	
	bool plot_write_b9dec1ef4eda5b08a0137ef1566a1fbd(class ::stat_tool::Convolution const & self_795c302c81665a5e9b66593bb34d8618, class ::stat_tool::StatError & error, std::string const & prefix, std::string const & title);
	
	bool plot_print_eb16d8e35f49592fb427e566f41924a5(class ::stat_tool::Distribution const & self_ba297ade4ec2559b8c17e36ec5555306, std::string const & path, class ::stat_tool::FrequencyDistribution const * histo);
	
	bool ascii_write_f398cc61e94d568786edf2866ab347c6(class ::stat_tool::DiscreteMixture const & self_9839f8684cb95ec091ee7dde51a7e147, class ::stat_tool::StatError & error, std::string const & path, bool exhaustive);
	
	bool plot_write_a992ac2d1c5a59cbb78dcd8f3183b1ea(class ::stat_tool::Compound const & self_e95f705e6c0c58008e5a01fae73debe6, class ::stat_tool::StatError & error, std::string const & prefix, std::string const & title);
	
	bool spreadsheet_write_8ad9deeff39f520b8db0f228142d0b86(class ::stat_tool::DistanceMatrix const & self_44dac8db1ff45337b3c1be9d33e36b56, class ::stat_tool::StatError & error, std::string const & path);
	
	bool spreadsheet_write_072f1eb9e79b5ec4ac0f0504c4bd2d98(class ::stat_tool::CompoundData const & self_3df29b17338354d49cd8f822201e189f, class ::stat_tool::StatError & error, std::string const & path);
	
	bool spreadsheet_write_90946d95145752f2875d667345f920c5(class ::stat_tool::DiscreteDistributionData const & self_699cdb7e19d058c5bb8d30066fdc55ac, class ::stat_tool::StatError & error, std::string const & path);
	
	bool plot_write_d41792e7c68c5f7a9fe968a4087e9534(class ::stat_tool::DiscreteMixture const & self_9839f8684cb95ec091ee7dde51a7e147, class ::stat_tool::StatError & error, std::string const & prefix, std::string const & title);
	
	bool ascii_write_61ed9e47ffda5e4ca7aae6f3938c8744(class ::stat_tool::MixtureData const & self_0bd75707ac98530d926ac17eb4fe47d3, class ::stat_tool::StatError & error, std::string const & path, bool exhaustive);
	
	std::string get_label_565170fe479d5248a4b31c8011783af4(class ::stat_tool::DistanceMatrix const & self_44dac8db1ff45337b3c1be9d33e36b56);
	
	bool plot_print_a13bbc6b2c2357d498d9e5b06d85bd3c(class ::stat_tool::ContinuousParametric & self_16995d5158735f999b84ef5efdd8439b, std::string const & path, class ::stat_tool::Histogram const * histo1, class ::stat_tool::FrequencyDistribution const * histo2);
	
	class ::stat_tool::DiscreteDistributionData * ascii_read_d00b270bd3835f568531951506bb4410(class ::stat_tool::StatError & error, std::string const & path);
	
	class ::stat_tool::MultivariateMixture * ascii_read_dd294317ade456caa8dd28e30146bae0(class ::stat_tool::StatError & error, std::string const & path, double cumul_threshold);
	
	bool spreadsheet_write_8553dbeb25f9567caf5a811f2bcf142e(class ::stat_tool::ConvolutionData const & self_acd6178a21305dc0a048b18ff459a4b3, class ::stat_tool::StatError & error, std::string const & path);
	
	class ::stat_tool::DiscreteParametricModel * ascii_read_d03ab63fd3665f35b2297d3df01fe531(class ::stat_tool::StatError & error, std::string const & path, double cumul_threshold);
	
	bool ascii_data_write_6a677fd70c1251e1bee86ad7330b2352(class ::stat_tool::Vectors const & self_e00cc711eccf595f84c6952a456d7688, class ::stat_tool::StatError & error, std::string const & path, bool exhaustive);
	
	bool plot_print_df2efb699a3a5b4486c73e45d656e875(class ::stat_tool::Curves const & self_42d9fc9d560e519a8cee703c630b1d4b, std::string const & path, int ilength, class ::stat_tool::Curves const * curves_0, class ::stat_tool::Curves const * curves_1);
	
	bool contingency_table_53063479d90a5bb5bad885693e2e5f52(class ::stat_tool::Vectors const & self_e00cc711eccf595f84c6952a456d7688, class ::stat_tool::StatError & error, class ::std::basic_ostream<char, std::char_traits<char> > & os, int variable1, int variable2, std::string const & path, enum ::stat_tool::output_format format);
	
	bool plot_write_5258db446d425dfa946d1bc7102f5ad6(class ::stat_tool::StatInterface const & self_57e8416c0462569e8916250db94d4d5a, class ::stat_tool::StatError & error, std::string const & prefix, std::string const & title);
	
	bool plot_print_46100a53a5955a16a690a06446f7bef8(class ::stat_tool::FrequencyDistribution const & self_f7aa7e5718795bbe854edb6369819dd9, std::string const & path, int nb_histo, class ::stat_tool::FrequencyDistribution const * * histo);
	
	bool selected_identifier_checking_245c9d680198502aba055a6484495a91(class ::stat_tool::StatError & error, int nb_individual, int * identifier, int nb_selected_individual, int * selected_identifier, std::string const & data_label);
	
	bool plot_write_b189acfb352a56228e3381403ac7c124(class ::stat_tool::VectorDistance const & self_01172f25879c56bd9d3f985ccab12de5, class ::stat_tool::StatError & error, std::string const & prefix, std::string const & title);
	
	class ::stat_tool::VectorDistance * ascii_read_9af48b1b53c3587ebda950257541f71d(class ::stat_tool::StatError & error, std::string const & path);
	
	bool spreadsheet_write_606ac800ccf550858b6d60ec81c6e17e(class ::stat_tool::Dendrogram const & self_acbc6c0efe495dce91ea3a9523ac2fa9, class ::stat_tool::StatError & error, std::string const & path);
	
	bool plot_write_faea833d86635b5ab2ba8dd2ecb423e5(class ::stat_tool::Distribution const & self_ba297ade4ec2559b8c17e36ec5555306, class ::stat_tool::StatError & error, std::string const & prefix, int nb_dist, class ::stat_tool::Distribution const * * idist, std::string const & title);
	
	bool ascii_write_becffa669e6e5eda9f9927803f6ab201(class ::stat_tool::DiscreteParametricModel const & self_96585eb2e8955e26852987df1de3cbd1, class ::stat_tool::StatError & error, std::string const & path, bool exhaustive);
	
	bool spreadsheet_write_27d18893029750f89f47eb35cf97486f(class ::stat_tool::DiscreteParametricModel const & self_96585eb2e8955e26852987df1de3cbd1, class ::stat_tool::StatError & error, std::string const & path);
	
	bool survival_plot_write_51f06af805d55e36a01c2d72302a95ac(class ::stat_tool::Distribution const & self_ba297ade4ec2559b8c17e36ec5555306, class ::stat_tool::StatError & error, std::string const & prefix, std::string const & title);
	
	class ::stat_tool::Vectors * ascii_read_53b245b9a04c5c53a163bd46d573d643(class ::stat_tool::StatError & error, std::string const & path);
	
	bool plot_print_78d1a7a202bf5eadbd531ee4c01557da(class ::stat_tool::Histogram const & self_73bc1367ebc25e6da2e7db8af7e3d828, std::string const & path);
	
	bool ascii_write_0d934e54b58f523aac0af7ab3df12737(class ::stat_tool::DistanceMatrix const & self_44dac8db1ff45337b3c1be9d33e36b56, class ::stat_tool::StatError & error, std::string const & path, bool exhaustive);
	
	bool ascii_write_814d58feea67507aa916b48b69053849(class ::stat_tool::Vectors const & self_e00cc711eccf595f84c6952a456d7688, class ::stat_tool::StatError & error, std::string const & path, bool exhaustive);
	
	bool ascii_data_write_7ab955beee015ce3aa4deb68aa6253a6(class ::stat_tool::MixtureData const & self_0bd75707ac98530d926ac17eb4fe47d3, class ::stat_tool::StatError & error, std::string const & path, bool exhaustive);
	
	bool survival_ascii_write_c4f4150db950536390b85cdd7c2439d1(class ::stat_tool::FrequencyDistribution const & self_f7aa7e5718795bbe854edb6369819dd9, class ::stat_tool::StatError & error, std::string const & path);
	
	bool spreadsheet_write_1c5a66445139510c98fbea38668616c3(class ::stat_tool::Regression const & self_6d51b373401f53d4aa43280d50ff7e29, class ::stat_tool::StatError & error, std::string const & path);
	
	bool plot_write_3f9868fc6a74577ba4375ee5c0513f63(class ::stat_tool::DiscreteMixtureData const & self_56a23f2d57545fffaaa0077ef13ef0ff, class ::stat_tool::StatError & error, std::string const & prefix, std::string const & title);
	
	bool plot_write_736abca042df51ee91ebcd6bb8d1352b(class ::stat_tool::Vectors const & self_e00cc711eccf595f84c6952a456d7688, class ::stat_tool::StatError & error, std::string const & prefix, std::string const & title);
	
	bool ascii_write_b5190834e6025b0891440384ac10d79f(class ::stat_tool::ConvolutionData const & self_acd6178a21305dc0a048b18ff459a4b3, class ::stat_tool::StatError & error, std::string const & path, bool exhaustive);
	
	bool plot_write_389cc509970958b79e41acd4f7777d94(class ::stat_tool::FrequencyDistribution const & self_f7aa7e5718795bbe854edb6369819dd9, class ::stat_tool::StatError & error, std::string const & prefix, int nb_histo, class ::stat_tool::FrequencyDistribution const * * ihisto, std::string const & title);
	
	class ::stat_tool::Convolution * ascii_read_63eca60a00ea53b1b2d40b39b77ff132(class ::stat_tool::StatError & error, std::string const & path, double cumul_threshold);
	
	bool ascii_write_0e14e13e449353cc91ac76673b04c770(class ::stat_tool::FrequencyDistribution const & self_f7aa7e5718795bbe854edb6369819dd9, class ::stat_tool::StatError & error, std::string const & path);
	
	bool plot_print_67285702bbb45140b7f6f669c80bcb44(class ::stat_tool::CategoricalProcess const & self_de7886a6095e5fde840400e992e24385, std::string const & prefix, std::string const & title, int process, class ::stat_tool::FrequencyDistribution * * empirical_observation, class ::stat_tool::FrequencyDistribution * marginal_distribution, enum ::stat_tool::model_type model);
	
	bool plot_write_fa95367df20e5171b515d01e57c83a5f(class ::stat_tool::ConvolutionData const & self_acd6178a21305dc0a048b18ff459a4b3, class ::stat_tool::StatError & error, std::string const & prefix, std::string const & title);
	
	bool plot_write_0e2dfd785bc55ed1b1860da0da6dfed2(class ::stat_tool::DistanceMatrix const & self_44dac8db1ff45337b3c1be9d33e36b56, class ::stat_tool::StatError & error, std::string const & prefix, std::string const & title);
	
	bool survival_spreadsheet_write_2d5c6910e527599dbd0d59aaca914acb(class ::stat_tool::Distribution const & self_ba297ade4ec2559b8c17e36ec5555306, class ::stat_tool::StatError & error, std::string const & path);
	
	bool plot_write_a9ee3ffd18fe5e4bbf35291c48ec6e3a(class ::stat_tool::DiscreteParametricModel const & self_96585eb2e8955e26852987df1de3cbd1, class ::stat_tool::StatError & error, std::string const & prefix, std::string const & title);
	
	class ::stat_tool::Compound * ascii_read_9c14f3a87c745f92af4f91bc1630e570(class ::stat_tool::StatError & error, std::string const & path, double cumul_threshold);
	
	bool plot_write_f2d5fcb84fbd5e5c8a7c72ade72be7a8(class ::stat_tool::Dendrogram const & self_acbc6c0efe495dce91ea3a9523ac2fa9, class ::stat_tool::StatError & error, std::string const & prefix, std::string const & title);
	
	bool plot_write_5e0c65eeae045977a3bb1aa32ed7e454(class ::stat_tool::MultivariateMixture const & self_2039c7a27ce4528da8cfc9e7ca6c4789, class ::stat_tool::StatError & error, std::string const & prefix, std::string const & title);
	
	void update_d75705f9274e5f8fb58207470a2fe6f4(class ::stat_tool::StatError & self_4cdf3268c85b531ab98dd99593a41730, std::string const & ilabel, int iline, int icolumn);
	
	bool hierarchical_clustering_9a8e184bcb2854c49e5153eabb514713(class ::stat_tool::DistanceMatrix const & self_44dac8db1ff45337b3c1be9d33e36b56, class ::stat_tool::StatError & error, class ::std::basic_ostream<char, std::char_traits<char> > & os, enum ::stat_tool::hierarchical_strategy strategy, enum ::stat_tool::linkage criterion, std::string const & path, enum ::stat_tool::output_format format);
	
	bool plot_write_630fefd35a7258fbbcac869efee69426(class ::stat_tool::CompoundData const & self_3df29b17338354d49cd8f822201e189f, class ::stat_tool::StatError & error, std::string const & prefix, std::string const & title);
	
	bool spreadsheet_write_b185ee74bac3528cb8dedace5dbe0c90(class ::stat_tool::VectorDistance const & self_01172f25879c56bd9d3f985ccab12de5, class ::stat_tool::StatError & error, std::string const & path);
	
	bool plot_print_3ac1e584c77a5bfdb16843bbe1851fed(class ::stat_tool::RegressionKernel const & self_ee3ef581666a5a45ab4576cfe45beb11, std::string const & path);
	
	bool spreadsheet_write_59596dfca7925fd89a6d320b8321a712(class ::stat_tool::MultivariateMixtureData const & self_5cde05a9631d5dcebc04732a82e9a74e, class ::stat_tool::StatError & error, std::string const & path);
	
	bool spreadsheet_write_2205836a245353dfbc9950f8b6f41628(class ::stat_tool::Vectors const & self_e00cc711eccf595f84c6952a456d7688, class ::stat_tool::StatError & error, std::string const & path);
	
	bool ascii_write_85d2ffa58a785b23be0e7992162c00c6(class ::stat_tool::DiscreteMixtureData const & self_56a23f2d57545fffaaa0077ef13ef0ff, class ::stat_tool::StatError & error, std::string const & path, bool exhaustive);
	
	bool spreadsheet_write_e8dd18fc792e54e1b0bd21e45c8951cb(class ::stat_tool::Mixture const & self_7fd12c45b7445cbda440f447ba16a244, class ::stat_tool::StatError & error, std::string const & path);
	
	bool ascii_write_898d9b4a1eb155b2b0d3e6fd81d0d5cd(class ::stat_tool::CompoundData const & self_3df29b17338354d49cd8f822201e189f, class ::stat_tool::StatError & error, std::string const & path, bool exhaustive);
	
	bool survival_plot_print_3e1ac3f9eb395e2bb3bc73d6e951e3aa(class ::stat_tool::FrequencyDistribution const & self_f7aa7e5718795bbe854edb6369819dd9, std::string const & path, double * survivor);
	
	bool plot_print_f332c698fc295076af58056713ddc526(std::string const & path, int nb_dist, class ::stat_tool::Distribution const * * dist, double * scale, int * dist_nb_value, int nb_histo, class ::stat_tool::FrequencyDistribution const * * histo);
	
	bool survival_plot_write_ed7d4d6631ea5a8d80985146babbde33(class ::stat_tool::FrequencyDistribution const & self_f7aa7e5718795bbe854edb6369819dd9, class ::stat_tool::StatError & error, std::string const & prefix, std::string const & title);
	
	bool plot_print_standard_residual_f58b54bba7b05235878222d98ddc2633(class ::stat_tool::Curves const & self_42d9fc9d560e519a8cee703c630b1d4b, std::string const & path, double * standard_residual);
	
	bool ascii_write_96a156ac7dd757d2926081b833167df1(class ::stat_tool::Regression const & self_6d51b373401f53d4aa43280d50ff7e29, class ::stat_tool::StatError & error, std::string const & path, bool exhaustive);
	
	class ::stat_tool::DiscreteMixture * ascii_read_543f048128b55b8ca39ac22bf2f802ce(class ::stat_tool::StatError & error, std::string const & path, double cumul_threshold);
	
	bool plot_write_6357a3ce5598512baa4fa593427a4fe2(class ::stat_tool::Clusters const & self_cf43b79f5a78554dace960313eb9c09b, class ::stat_tool::StatError & error, std::string const & prefix, std::string const & title);
	
	bool ascii_write_8ac6af2e6b0b52ee9c1cc81342a8a8d6(class ::stat_tool::DiscreteDistributionData const & self_699cdb7e19d058c5bb8d30066fdc55ac, class ::stat_tool::StatError & error, std::string const & path, bool exhaustive);
	
	bool rank_correlation_computation_00d124e24d345c9f927ea5d72a2a405e(class ::stat_tool::Vectors const & self_e00cc711eccf595f84c6952a456d7688, class ::stat_tool::StatError & error, class ::std::basic_ostream<char, std::char_traits<char> > & os, enum ::stat_tool::correlation_type correl_type, std::string const & path);
	
	bool dissimilarity_ascii_write_3ed836e2b03f5396a9daafcc5fa9e96f(class ::stat_tool::FrequencyDistribution const & self_f7aa7e5718795bbe854edb6369819dd9, class ::stat_tool::StatError & error, std::string const & path, int nb_histo, class ::stat_tool::FrequencyDistribution const * * ihisto, enum ::stat_tool::variable_type type, double * * dissimilarity);
	
	bool ascii_write_4e98d150f1905ef5b4bb5fde8e2f7715(class ::stat_tool::Clusters const & self_cf43b79f5a78554dace960313eb9c09b, class ::stat_tool::StatError & error, std::string const & path, bool exhaustive);
	
	bool comparison_9f04deaf7bdb51f4804495d5a3f4c7dd(class ::stat_tool::FrequencyDistribution const & self_f7aa7e5718795bbe854edb6369819dd9, class ::stat_tool::StatError & error, class ::std::basic_ostream<char, std::char_traits<char> > & os, int nb_histo, class ::stat_tool::FrequencyDistribution const * * ihisto, enum ::stat_tool::variable_type type, std::string const & path, enum ::stat_tool::output_format format);
	
	void correction_update_5ac99783e68f5be88e3c29ca5f6cf0dd(class ::stat_tool::StatError & self_4cdf3268c85b531ab98dd99593a41730, std::string const & ilabel, int correction, int iline, int icolumn);
	
	bool dissimilarity_spreadsheet_write_149cbcd3090956fdbd4ea073b8e746af(class ::stat_tool::FrequencyDistribution const & self_f7aa7e5718795bbe854edb6369819dd9, class ::stat_tool::StatError & error, std::string const & path, int nb_histo, class ::stat_tool::FrequencyDistribution const * * ihisto, enum ::stat_tool::variable_type type, double * * dissimilarity);
	
	bool spreadsheet_write_f67b5bb07dd553a7b97be9ac9024e0d3(class ::stat_tool::DiscreteMixtureData const & self_56a23f2d57545fffaaa0077ef13ef0ff, class ::stat_tool::StatError & error, std::string const & path);
	
	bool plot_write_40e4d8a6d6e55643a6ba930122dba403(class ::stat_tool::MultivariateMixtureData const & self_5cde05a9631d5dcebc04732a82e9a74e, class ::stat_tool::StatError & error, std::string const & prefix, std::string const & title);
	
	bool plot_print_23e763005a63572c8caa6b80b38134f0(class ::stat_tool::ContinuousParametricProcess const & self_9ba34b38b8a75c968e25c23f75118106, std::string const & prefix, std::string const & title, int process, class ::stat_tool::Histogram * * observation_histogram, class ::stat_tool::FrequencyDistribution * * observation_distribution, class ::stat_tool::Histogram * marginal_histogram, class ::stat_tool::FrequencyDistribution * marginal_distribution, int nb_value, double * * empirical_cdf, enum ::stat_tool::model_type model);
	
	bool ascii_write_ed6534b1b5cb502dbcd77e2db66484d7(class ::stat_tool::Mixture const & self_7fd12c45b7445cbda440f447ba16a244, class ::stat_tool::StatError & error, std::string const & path, bool exhaustive);
	
	bool plot_write_847e51c8701752b8940f2ed779319379(class ::stat_tool::Regression const & self_6d51b373401f53d4aa43280d50ff7e29, class ::stat_tool::StatError & error, std::string const & prefix, std::string const & title);
	
	bool plot_print_f302824b98df5aebba39caa19d1535e0(class ::stat_tool::Distribution const & self_ba297ade4ec2559b8c17e36ec5555306, std::string const & path, double * concentration, double scale);
	
	bool survival_ascii_write_e88be3ac861954a7945656be4af3766c(class ::stat_tool::Distribution const & self_ba297ade4ec2559b8c17e36ec5555306, class ::stat_tool::StatError & error, std::string const & path);
	
	bool plot_print_267ccb22b12b5467b535ef445643b047(class ::stat_tool::DiscreteParametricProcess const & self_003de01bc70a5a99867c38bed57c58a5, std::string const & prefix, std::string const & title, int process, class ::stat_tool::FrequencyDistribution * * empirical_observation, class ::stat_tool::FrequencyDistribution * marginal_distribution, enum ::stat_tool::model_type model);
	
	bool ascii_write_f5455e0a03ce56edb9affad684476e91(class ::stat_tool::VectorDistance const & self_01172f25879c56bd9d3f985ccab12de5, class ::stat_tool::StatError & error, std::string const & path, bool exhaustive);
	
	bool ascii_write_b416f25d3ac95ab1bb8a73e51d624831(class ::stat_tool::MultivariateMixture const & self_2039c7a27ce4528da8cfc9e7ca6c4789, class ::stat_tool::StatError & error, std::string const & path, bool exhaustive);
	
	bool q_q_plot_print_635196db3f585030b33fb62ba232e490(class ::stat_tool::ContinuousParametric const & self_16995d5158735f999b84ef5efdd8439b, std::string const & path, int nb_value, double * * empirical_cdf);
	
	bool plot_write_2854adbc735d583483d3fff059d40803(class ::stat_tool::MixtureData const & self_0bd75707ac98530d926ac17eb4fe47d3, class ::stat_tool::StatError & error, std::string const & prefix, std::string const & title);
	
	bool ascii_write_999c19aa91ed56cd8f4dd0fb6d325e48(class ::stat_tool::Dendrogram const & self_acbc6c0efe495dce91ea3a9523ac2fa9, class ::stat_tool::StatError & error, std::string const & path, bool exhaustive);
	
	bool ascii_write_d94a583fc98b5c01baba7a2c77c562a9(class ::stat_tool::Convolution const & self_795c302c81665a5e9b66593bb34d8618, class ::stat_tool::StatError & error, std::string const & path, bool exhaustive);
	
	class ::stat_tool::Mixture * ascii_read_fcac4ba71c2b519e85e7fcba6c2596a3(class ::stat_tool::StatError & error, std::string const & path, double cumul_threshold);
	
	bool spreadsheet_write_a612b81f8c265e589b9dc2c24a927fe6(class ::stat_tool::Clusters const & self_cf43b79f5a78554dace960313eb9c09b, class ::stat_tool::StatError & error, std::string const & path);
	
	bool spreadsheet_write_818d7344e8755d2bb2f8f377a8c61d3d(class ::stat_tool::MultivariateMixture const & self_2039c7a27ce4528da8cfc9e7ca6c4789, class ::stat_tool::StatError & error, std::string const & path);
	
	bool spreadsheet_write_b486a5f6069e5ee1b978538a2f8fee74(class ::stat_tool::MixtureData const & self_0bd75707ac98530d926ac17eb4fe47d3, class ::stat_tool::StatError & error, std::string const & path);
	
	bool spreadsheet_write_8c5a70e9b2bd50b090deec3cb3a905c3(class ::stat_tool::Compound const & self_e95f705e6c0c58008e5a01fae73debe6, class ::stat_tool::StatError & error, std::string const & path);
	
	bool variance_analysis_83315e4f4ac958a9afa637694b780c11(class ::stat_tool::Vectors const & self_e00cc711eccf595f84c6952a456d7688, class ::stat_tool::StatError & error, class ::std::basic_ostream<char, std::char_traits<char> > & os, int class_variable, int response_variable, int response_type, std::string const & path, enum ::stat_tool::output_format format);
	
	bool survival_spreadsheet_write_e9fb2b43f944590a9adc613eeece2e69(class ::stat_tool::FrequencyDistribution const & self_f7aa7e5718795bbe854edb6369819dd9, class ::stat_tool::StatError & error, std::string const & path);
	
	bool plot_print_2200fedbc57a5e2a96d155e65fd63cb6(class ::stat_tool::FrequencyDistribution const & self_f7aa7e5718795bbe854edb6369819dd9, std::string const & path, double * cumul, double * concentration, double shift);
	
	bool survival_plot_print_4817f9e42ffc5c72bf8b7dea51d5a044(class ::stat_tool::Distribution const & self_ba297ade4ec2559b8c17e36ec5555306, std::string const & path, double * survivor);
	
	bool ascii_write_586d21963d015068bdf7f31d8d39479c(class ::stat_tool::StatInterface const & self_57e8416c0462569e8916250db94d4d5a, class ::stat_tool::StatError & error, std::string const & path, bool exhaustive);
	
	bool spreadsheet_write_4b9baa0039c75c6aaa2d3c2a11210cb9(class ::stat_tool::DiscreteMixture const & self_9839f8684cb95ec091ee7dde51a7e147, class ::stat_tool::StatError & error, std::string const & path);
	
	bool spreadsheet_write_e5966fbb35385a33a0448c8a567a2c02(class ::stat_tool::Convolution const & self_795c302c81665a5e9b66593bb34d8618, class ::stat_tool::StatError & error, std::string const & path);
	
	bool ascii_write_9a47c2d60b6b5455b0bc51489402612a(class ::stat_tool::Compound const & self_e95f705e6c0c58008e5a01fae73debe6, class ::stat_tool::StatError & error, std::string const & path, bool exhaustive);
}