#if ( defined(DIPOLS) || defined(ONE_CELL_SN) || defined(INSERTSNONCEASTEP) ) && !defined(SNE_DISTR)
#error SNE_DISTR must be defined in this configuration!
#endif

#if ( defined(MASS_COMPENS) || defined(VZ_LIMITS) ) && !defined(ANY_LIMITS)
#error ANY_LIMITS must be defined in this configuration!
#endif

#if defined(DIPOLS) && !defined(MAGNETIC)
#error MAGNETIC must be defined while defined DIPOLS
#endif

#if defined(SN_DISTRIBUTION) && !defined(SNE_DISTR)
#error SNE_DISTR must be defined while defined SN_DISTRIBUTION
#endif

#if defined(DISTR_INFLOW) && !defined(MASS_COMPENS)
#error MASS_COMPENS must be defined while defined DISTR_INFLOW
#endif
