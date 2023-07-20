/*! \file config.h
 *
 *  \brief Contains macro definitions.
 *         Also reserved for future `configure` script
 *
 */
#ifndef VERBOSE_LEVEL
#  define VERBOSE_LEVEL 1
#endif

#if VERBOSE_LEVEL > 3
#  define INFO_ENTRY(fname) print*,"[<] ",fname," (): entering"
#  define INFO_EXIT(fname)  print*,"[>] ",fname," (): exiting"
#else
#  define INFO_ENTRY(fname) 
#  define INFO_EXIT(fname) 
#endif

/* The used Git version */
#define XSTR(VERS) #VERS
#define STR(VERS) XSTR(VERS)

/* possible values of the evolution_mode */
#define EM_NSE     1
#define EM_NETHOT  2
#define EM_NETCOLD 3
#define EM_DECAYS  4
#define EM_INIT    5
#define EM_TERMINATE 6

/* possible values for the space of the interpolation */
#define LINLIN     1
#define LINLOG     2
#define LOGLIN     3
#define LOGLOG     4


/* possible values for the type of the interpolation */
#define itype_LINEAR     1
#define itype_CUBIC      2
#define itype_AKIMA      3
#define itype_MAKIMA     4
#define itype_PCHIP      5

/* define real kind to doubles */
#define r_kind     8

/* possible values of reaction sources                    */
#define rrs_reacl 1 /* Reaclib reactions                  */
#define rrs_tabl  2 /* Tabulated reactions                */
#define rrs_nu    3 /* Neutrino reactions                 */
#define rrs_twr   4 /* Theoretical weak rate              */
#define rrs_detb  5 /* Calculated detailed balance rate   */
#define rrs_fiss  6 /* Fission reaction                   */
#define rrs_wext  7 /* Beta delayed neutron emission rate */
#define rrs_aext  8 /* Additional alpha decay rate        */

/* possible values of reaction rate typ */
#define rrt_betm   1 /* Beta-minus                */
#define rrt_betp   2 /* Beta-plus                 */
#define rrt_ec     3 /* Electron capture          */
#define rrt_alpd   4 /* Alpha-decay               */
#define rrt_pemi   5 /* Proton-emission           */
#define rrt_nemi   6 /* Neutron-emission          */
#define rrt_ng     7 /* (n,gamma)                 */
#define rrt_gn     8 /* (gamma,n)                 */
#define rrt_ag     9 /* (alpha,gamma)             */
#define rrt_ga    10 /* (gamma,alpha)             */
#define rrt_pg    11 /* (p,gamma)                 */
#define rrt_gp    12 /* (gamma,p)                 */
#define rrt_na    13 /* (n,alpha)                 */
#define rrt_an    14 /* (alpha,n)                 */
#define rrt_np    15 /* (n,p)                     */
#define rrt_pn    16 /* (p,n)                     */
#define rrt_pa    17 /* (p,alpha)                 */
#define rrt_ap    18 /* (alpha,p)                 */
#define rrt_nu    19 /* Neutrino                  */
#define rrt_nf    20 /* neutron induced fission   */
#define rrt_bf    21 /* beta-delayed fission      */
#define rrt_sf    22 /* spontaneous fission       */
#define rrt_o     23 /* other                     */


