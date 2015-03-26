#ifndef __CLAPACK_H
#define __CLAPACK_H

#include "numerics/f2c.h"

#ifdef NO_UNDERSCORE
#define F77_CALL(x) x
#else
#define F77_CALL(x) x ## _
#endif


#ifdef F77_UPPER
#define dbdsdc DBDSDC
#define dbdsqr DBDSQR
#define ddisna DDISNA
#define dgbbrd DGBBRD
#define dgbcon DGBCON
#define dgbequ DGBEQU
#define dgbrfs DGBRFS
#define dgbsv DGBSV
#define dgbsvx DGBSVX
#define dgbtrf DGBTRF
#define dgbtrs DGBTRS
#define dgebak DGEBAK
#define dgebal DGEBAL
#define dgebrd DGEBRD
#define dgecon DGECON
#define dgeequ DGEEQU
#define dgees DGEES
#define dgeesx DGEESX
#define dgeev DGEEV
#define dgeevx DGEEVX
#define dgegs DGEGS
#define dgegv DGEGV
#define dgehrd DGEHRD
#define dgelqf DGELQF
#define dgels DGELS
#define dgelsd DGELSD
#define dgelss DGELSS
#define dgelsx DGELSX
#define dgelsy DGELSY
#define dgeqlf DGEQLF
#define dgeqpf DGEQPF
#define dgeqrf DGEQRF
#define dgerfs DGERFS
#define dgerqf DGERQF
#define dgesdd DGESDD
#define dgesv DGESV
#define dgesvd DGESVD
#define dgesvx DGESVX
#define dgetrf DGETRF
#define dgetri DGETRI
#define dgetrs DGETRS
#define dggbak DGGBAK
#define dggbal DGGBAL
#define dgges DGGES
#define dggesx DGGESX
#define dggev DGGEV
#define dggevx DGGEVX
#define dggglm DGGGLM
#define dgghrd DGGHRD
#define dgglse DGGLSE
#define dggqrf DGGQRF
#define dggrqf DGGRQF
#define dggsvd DGGSVD
#define dggsvp DGGSVP
#define dgtcon DGTCON
#define dgtrfs DGTRFS
#define dgtsv DGTSV
#define dgtsvx DGTSVX
#define dgttrf DGTTRF
#define dgttrs DGTTRS
#define dhgeqz DHGEQZ
#define dhsein DHSEIN
#define dhseqr DHSEQR
#define dlabad DLABAD
#define dlabrd DLABRD
#define dlacon DLACON
#define dlacpy DLACPY
#define dladiv DLADIV
#define dlaebz DLAEBZ
#define dlaeda DLAEDA
#define dlaein DLAEIN
#define dlaexc DLAEXC
#define dlagtf DLAGTF
#define dlagtm DLAGTM
#define dlagts DLAGTS
#define dlahqr DLAHQR
#define dlahrd DLAHRD
#define dlalsa DLALSA
#define dlalsd DLALSD
#define dlamrg DLAMRG
#define dlapll DLAPLL
#define dlapmt DLAPMT
#define dlaqgb DLAQGB
#define dlaqge DLAQGE
#define dlaqps DLAQPS
#define dlaqsb DLAQSB
#define dlaqsp DLAQSP
#define dlaqsy DLAQSY
#define dlaqtr DLAQTR
#define dlarf DLARF
#define dlarfb DLARFB
#define dlarfg DLARFG
#define dlarft DLARFT
#define dlarfx DLARFX
#define dlargv DLARGV
#define dlarnv DLARNV
#define dlarrb DLARRB
#define dlarre DLARRE
#define dlarrf DLARRF
#define dlarrv DLARRV
#define dlartg DLARTG
#define dlartv DLARTV
#define dlaruv DLARUV
#define dlarz DLARZ
#define dlarzb DLARZB
#define dlarzt DLARZT
#define dlascl DLASCL
#define dlasda DLASDA
#define dlasdq DLASDQ
#define dlasdt DLASDT
#define dlaset DLASET
#define dlasr DLASR
#define dlasrt DLASRT
#define dlassq DLASSQ
#define dlaswp DLASWP
#define dlasyf DLASYF
#define dlatbs DLATBS
#define dlatdf DLATDF
#define dlatps DLATPS
#define dlatrd DLATRD
#define dlatrs DLATRS
#define dlatrz DLATRZ
#define dlatzm DLATZM
#define dlauum DLAUUM
#define dopgtr DOPGTR
#define dopmtr DOPMTR
#define dorgbr DORGBR
#define dorghr DORGHR
#define dorglq DORGLQ
#define dorgql DORGQL
#define dorgqr DORGQR
#define dorgrq DORGRQ
#define dorgtr DORGTR
#define dormbr DORMBR
#define dormhr DORMHR
#define dormlq DORMLQ
#define dormql DORMQL
#define dormqr DORMQR
#define dormrq DORMRQ
#define dormrz DORMRZ
#define dormtr DORMTR
#define dpbcon DPBCON
#define dpbequ DPBEQU
#define dpbrfs DPBRFS
#define dpbstf DPBSTF
#define dpbsv DPBSV
#define dpbsvx DPBSVX
#define dpbtrf DPBTRF
#define dpbtrs DPBTRS
#define dpocon DPOCON
#define dpoequ DPOEQU
#define dporfs DPORFS
#define dposv DPOSV
#define dposvx DPOSVX
#define dpotrf DPOTRF
#define dpotri DPOTRI
#define dpotrs DPOTRS
#define dppcon DPPCON
#define dppequ DPPEQU
#define dpprfs DPPRFS
#define dppsv DPPSV
#define dppsvx DPPSVX
#define dpptrf DPPTRF
#define dpptri DPPTRI
#define dpptrs DPPTRS
#define dptcon DPTCON
#define dpteqr DPTEQR
#define dptrfs DPTRFS
#define dptsv DPTSV
#define dptsvx DPTSVX
#define dpttrf DPTTRF
#define dpttrs DPTTRS
#define drscl DRSCL
#define dsbev DSBEV
#define dsbevd DSBEVD
#define dsbevx DSBEVX
#define dsbgst DSBGST
#define dsbgv DSBGV
#define dsbgvd DSBGVD
#define dsbgvx DSBGVX
#define dsbtrd DSBTRD
#define dspcon DSPCON
#define dspev DSPEV
#define dspevd DSPEVD
#define dspevx DSPEVX
#define dspgst DSPGST
#define dspgv DSPGV
#define dspgvd DSPGVD
#define dspgvx DSPGVX
#define dsprfs DSPRFS
#define dspsv DSPSV
#define dspsvx DSPSVX
#define dsptrd DSPTRD
#define dsptrf DSPTRF
#define dsptri DSPTRI
#define dsptrs DSPTRS
#define dstebz DSTEBZ
#define dstedc DSTEDC
#define dstegr DSTEGR
#define dstein DSTEIN
#define dsteqr DSTEQR
#define dsterf DSTERF
#define dstev DSTEV
#define dstevd DSTEVD
#define dstevr DSTEVR
#define dstevx DSTEVX
#define dsycon DSYCON
#define dsyev DSYEV
#define dsyevd DSYEVD
#define dsyevr DSYEVR
#define dsyevx DSYEVX
#define dsygst DSYGST
#define dsygv DSYGV
#define dsygvd DSYGVD
#define dsygvx DSYGVX
#define dsyrfs DSYRFS
#define dsysv DSYSV
#define dsysvx DSYSVX
#define dsytrd DSYTRD
#define dsytrf DSYTRF
#define dsytri DSYTRI
#define dsytrs DSYTRS
#define dtbcon DTBCON
#define dtbrfs DTBRFS
#define dtbtrs DTBTRS
#define dtgevc DTGEVC
#define dtgexc DTGEXC
#define dtgsen DTGSEN
#define dtgsja DTGSJA
#define dtgsna DTGSNA
#define dtgsyl DTGSYL
#define dtpcon DTPCON
#define dtprfs DTPRFS
#define dtptri DTPTRI
#define dtptrs DTPTRS
#define dtrcon DTRCON
#define dtrevc DTREVC
#define dtrexc DTREXC
#define dtrrfs DTRRFS
#define dtrsen DTRSEN
#define dtrsna DTRSNA
#define dtrsyl DTRSYL
#define dtrtri DTRTRI
#define dtrtrs DTRTRS
#define dtzrqf DTZRQF
#define dtzrzf DTZRZF
#define ilaenv ILAENV
#define slasdt SLASDT
#define xerbla XERBLA
#endif




extern "C"
{
int F77_CALL(dbdsdc)(char *uplo, char *compq, int *n, double * d__, double *e, double *u, int *ldu, double *vt,  int *ldvt, double *q, int *iq, double *work, int * iwork, int *info);
int F77_CALL(dbdsqr)(char *uplo, int *n, int *ncvt, int * nru, int *ncc, double *d__, double *e, double *vt,  int *ldvt, double *u, int *ldu, double *c__, int * ldc, double *work, int *info);
int F77_CALL(ddisna)(char *job, int *m, int *n, double * d__, double *sep, int *info);
int F77_CALL(dgbbrd)(char *vect, int *m, int *n, int *ncc, int *kl, int *ku, double *ab, int *ldab, double * d__, double *e, double *q, int *ldq, double *pt,  int *ldpt, double *c__, int *ldc, double *work,  int *info);
int F77_CALL(dgbcon)(char *norm, int *n, int *kl, int *ku, double *ab, int *ldab, int *ipiv, double *anorm,  double *rcond, double *work, int *iwork, int *info);
int F77_CALL(dgbequ)(int *m, int *n, int *kl, int *ku, double *ab, int *ldab, double *r__, double *c__,  double *rowcnd, double *colcnd, double *amax, int * info);
int F77_CALL(dgbrfs)(char *trans, int *n, int *kl, int * ku, int *nrhs, double *ab, int *ldab, double *afb,  int *ldafb, int *ipiv, double *b, int *ldb,  double *x, int *ldx, double *ferr, double *berr,  double *work, int *iwork, int *info);
int F77_CALL(dgbsv)(int *n, int *kl, int *ku, int * nrhs, double *ab, int *ldab, int *ipiv, double *b,  int *ldb, int *info);
int F77_CALL(dgbsvx)(char *fact, char *trans, int *n, int *kl, int *ku, int *nrhs, double *ab, int *ldab,  double *afb, int *ldafb, int *ipiv, char *equed,  double *r__, double *c__, double *b, int *ldb,  double *x, int *ldx, double *rcond, double *ferr,  double *berr, double *work, int *iwork, int *info);
int F77_CALL(dgbtf2)(int *m, int *n, int *kl, int *ku, double *ab, int *ldab, int *ipiv, int *info);
int F77_CALL(dgbtrf)(int *m, int *n, int *kl, int *ku, double *ab, int *ldab, int *ipiv, int *info);
int F77_CALL(dgbtrs)(char *trans, int *n, int *kl, int * ku, int *nrhs, double *ab, int *ldab, int *ipiv,  double *b, int *ldb, int *info);
int F77_CALL(dgebak)(char *job, char *side, int *n, int *ilo,  int *ihi, double *scale, int *m, double *v, int * ldv, int *info);
int F77_CALL(dgebal)(char *job, int *n, double *a, int * lda, int *ilo, int *ihi, double *scale, int *info);
int F77_CALL(dgebd2)(int *m, int *n, double *a, int * lda, double *d__, double *e, double *tauq, double * taup, double *work, int *info);
int F77_CALL(dgebrd)(int *m, int *n, double *a, int * lda, double *d__, double *e, double *tauq, double * taup, double *work, int *lwork, int *info);
int F77_CALL(dgecon)(char *norm, int *n, double *a, int * lda, double *anorm, double *rcond, double *work, int * iwork, int *info);
int F77_CALL(dgeequ)(int *m, int *n, double *a, int * lda, double *r__, double *c__, double *rowcnd, double  *colcnd, double *amax, int *info);
int F77_CALL(dgees)(char *jobvs, char *sort, L_fp select, int *n,  double *a, int *lda, int *sdim, double *wr,  double *wi, double *vs, int *ldvs, double *work,  int *lwork, logical *bwork, int *info);
int F77_CALL(dgeesx)(char *jobvs, char *sort, L_fp select, char * sense, int *n, double *a, int *lda, int *sdim,  double *wr, double *wi, double *vs, int *ldvs,  double *rconde, double *rcondv, double *work, int * lwork, int *iwork, int *liwork, logical *bwork, int *info);
int F77_CALL(dgeev)(char *jobvl, char *jobvr, int *n, double * a, int *lda, double *wr, double *wi, double *vl,  int *ldvl, double *vr, int *ldvr, double *work,  int *lwork, int *info);
int F77_CALL(dgeevx)(char *balanc, char *jobvl, char *jobvr, char * sense, int *n, double *a, int *lda, double *wr,  double *wi, double *vl, int *ldvl, double *vr,  int *ldvr, int *ilo, int *ihi, double *scale,  double *abnrm, double *rconde, double *rcondv, double  *work, int *lwork, int *iwork, int *info);
int F77_CALL(dgegs)(char *jobvsl, char *jobvsr, int *n,  double *a, int *lda, double *b, int *ldb, double * alphar, double *alphai, double *beta, double *vsl,  int *ldvsl, double *vsr, int *ldvsr, double *work,  int *lwork, int *info);
int F77_CALL(dgegv)(char *jobvl, char *jobvr, int *n, double * a, int *lda, double *b, int *ldb, double *alphar,  double *alphai, double *beta, double *vl, int *ldvl,  double *vr, int *ldvr, double *work, int *lwork,  int *info);
int F77_CALL(dgehd2)(int *n, int *ilo, int *ihi,  double *a, int *lda, double *tau, double *work,  int *info);
int F77_CALL(dgehrd)(int *n, int *ilo, int *ihi,  double *a, int *lda, double *tau, double *work,  int *lwork, int *info);
int F77_CALL(dgelq2)(int *m, int *n, double *a, int * lda, double *tau, double *work, int *info);
int F77_CALL(dgelqf)(int *m, int *n, double *a, int * lda, double *tau, double *work, int *lwork, int *info);
int F77_CALL(dgels)(char *trans, int *m, int *n, int * nrhs, double *a, int *lda, double *b, int *ldb,  double *work, int *lwork, int *info);
int F77_CALL(dgelsd)(int *m, int *n, int *nrhs,  double *a, int *lda, double *b, int *ldb, double * s, double *rcond, int *rank, double *work, int *lwork, int *iwork, int *info);
int F77_CALL(dgelss)(int *m, int *n, int *nrhs,  double *a, int *lda, double *b, int *ldb, double * s, double *rcond, int *rank, double *work, int *lwork, int *info);
int F77_CALL(dgelsx)(int *m, int *n, int *nrhs,  double *a, int *lda, double *b, int *ldb, int * jpvt, double *rcond, int *rank, double *work, int * info);
int F77_CALL(dgelsy)(int *m, int *n, int *nrhs,  double *a, int *lda, double *b, int *ldb, int * jpvt, double *rcond, int *rank, double *work, int * lwork, int *info);
int F77_CALL(dgeql2)(int *m, int *n, double *a, int * lda, double *tau, double *work, int *info);
int F77_CALL(dgeqlf)(int *m, int *n, double *a, int * lda, double *tau, double *work, int *lwork, int *info);
int F77_CALL(dgeqp3)(int *m, int *n, double *a, int * lda, int *jpvt, double *tau, double *work, int *lwork, int *info);
int F77_CALL(dgeqpf)(int *m, int *n, double *a, int * lda, int *jpvt, double *tau, double *work, int *info);
int F77_CALL(dgeqr2)(int *m, int *n, double *a, int * lda, double *tau, double *work, int *info);
int F77_CALL(dgeqrf)(int *m, int *n, double *a, int * lda, double *tau, double *work, int *lwork, int *info);
int F77_CALL(dgerfs)(char *trans, int *n, int *nrhs,  double *a, int *lda, double *af, int *ldaf, int * ipiv, double *b, int *ldb, double *x, int *ldx,  double *ferr, double *berr, double *work, int *iwork,  int *info);
int F77_CALL(dgerq2)(int *m, int *n, double *a, int * lda, double *tau, double *work, int *info);
int F77_CALL(dgerqf)(int *m, int *n, double *a, int * lda, double *tau, double *work, int *lwork, int *info);
int F77_CALL(dgesc2)(int *n, double *a, int *lda,  double *rhs, int *ipiv, int *jpiv, double *scale);
int F77_CALL(dgesdd)(char *jobz, int *m, int *n, double * a, int *lda, double *s, double *u, int *ldu,  double *vt, int *ldvt, double *work, int *lwork,  int *iwork, int *info);
int F77_CALL(dgesv)(int *n, int *nrhs, double *a, int  *lda, int *ipiv, double *b, int *ldb, int *info);
int F77_CALL(dgesvd)(char *jobu, char *jobvt, int *m, int *n,  double *a, int *lda, double *s, double *u, int * ldu, double *vt, int *ldvt, double *work, int *lwork,  int *info);
int F77_CALL(dgesvx)(char *fact, char *trans, int *n, int * nrhs, double *a, int *lda, double *af, int *ldaf,  int *ipiv, char *equed, double *r__, double *c__,  double *b, int *ldb, double *x, int *ldx, double * rcond, double *ferr, double *berr, double *work, int * iwork, int *info);
int F77_CALL(dgetc2)(int *n, double *a, int *lda, int  *ipiv, int *jpiv, int *info);
int F77_CALL(dgetf2)(int *m, int *n, double *a, int * lda, int *ipiv, int *info);
int F77_CALL(dgetrf)(int *m, int *n, double *a, int * lda, int *ipiv, int *info);
int F77_CALL(dgetri)(int *n, double *a, int *lda, int  *ipiv, double *work, int *lwork, int *info);
int F77_CALL(dgetrs)(char *trans, int *n, int *nrhs,  double *a, int *lda, int *ipiv, double *b, int * ldb, int *info);
int F77_CALL(dggbak)(char *job, char *side, int *n, int *ilo,  int *ihi, double *lscale, double *rscale, int *m,  double *v, int *ldv, int *info);
int F77_CALL(dggbal)(char *job, int *n, double *a, int * lda, double *b, int *ldb, int *ilo, int *ihi,  double *lscale, double *rscale, double *work, int * info);
int F77_CALL(dgges)(char *jobvsl, char *jobvsr, char *sort, L_fp  delctg, int *n, double *a, int *lda, double *b,  int *ldb, int *sdim, double *alphar, double *alphai,  double *beta, double *vsl, int *ldvsl, double *vsr,  int *ldvsr, double *work, int *lwork, logical *bwork,  int *info);
int F77_CALL(dggesx)(char *jobvsl, char *jobvsr, char *sort, L_fp  delctg, char *sense, int *n, double *a, int *lda,  double *b, int *ldb, int *sdim, double *alphar,  double *alphai, double *beta, double *vsl, int *ldvsl, double *vsr, int *ldvsr, double *rconde, double * rcondv, double *work, int *lwork, int *iwork, int * liwork, logical *bwork, int *info);
int F77_CALL(dggev)(char *jobvl, char *jobvr, int *n, double * a, int *lda, double *b, int *ldb, double *alphar,  double *alphai, double *beta, double *vl, int *ldvl,  double *vr, int *ldvr, double *work, int *lwork,  int *info);
int F77_CALL(dggevx)(char *balanc, char *jobvl, char *jobvr, char * sense, int *n, double *a, int *lda, double *b,  int *ldb, double *alphar, double *alphai, double * beta, double *vl, int *ldvl, double *vr, int *ldvr,  int *ilo, int *ihi, double *lscale, double *rscale,  double *abnrm, double *bbnrm, double *rconde, double * rcondv, double *work, int *lwork, int *iwork, logical * bwork, int *info);
int F77_CALL(dggglm)(int *n, int *m, int *p, double * a, int *lda, double *b, int *ldb, double *d__,  double *x, double *y, double *work, int *lwork,  int *info);
int F77_CALL(dgghrd)(char *compq, char *compz, int *n, int * ilo, int *ihi, double *a, int *lda, double *b,  int *ldb, double *q, int *ldq, double *z__, int * ldz, int *info);
int F77_CALL(dgglse)(int *m, int *n, int *p, double * a, int *lda, double *b, int *ldb, double *c__,  double *d__, double *x, double *work, int *lwork,  int *info);
int F77_CALL(dggqrf)(int *n, int *m, int *p, double * a, int *lda, double *taua, double *b, int *ldb,  double *taub, double *work, int *lwork, int *info);
int F77_CALL(dggrqf)(int *m, int *p, int *n, double * a, int *lda, double *taua, double *b, int *ldb,  double *taub, double *work, int *lwork, int *info);
int F77_CALL(dggsvd)(char *jobu, char *jobv, char *jobq, int *m,  int *n, int *p, int *k, int *l, double *a,  int *lda, double *b, int *ldb, double *alpha,  double *beta, double *u, int *ldu, double *v, int  *ldv, double *q, int *ldq, double *work, int *iwork,  int *info);
int F77_CALL(dggsvp)(char *jobu, char *jobv, char *jobq, int *m,  int *p, int *n, double *a, int *lda, double *b,  int *ldb, double *tola, double *tolb, int *k, int  *l, double *u, int *ldu, double *v, int *ldv,  double *q, int *ldq, int *iwork, double *tau,  double *work, int *info);
int F77_CALL(dgtcon)(char *norm, int *n, double *dl,  double *d__, double *du, double *du2, int *ipiv,  double *anorm, double *rcond, double *work, int * iwork, int *info);
int F77_CALL(dgtrfs)(char *trans, int *n, int *nrhs,  double *dl, double *d__, double *du, double *dlf,  double *df, double *duf, double *du2, int *ipiv,  double *b, int *ldb, double *x, int *ldx, double * ferr, double *berr, double *work, int *iwork, int * info);
int F77_CALL(dgtsv)(int *n, int *nrhs, double *dl,  double *d__, double *du, double *b, int *ldb, int  *info);
int F77_CALL(dgtsvx)(char *fact, char *trans, int *n, int * nrhs, double *dl, double *d__, double *du, double * dlf, double *df, double *duf, double *du2, int *ipiv,  double *b, int *ldb, double *x, int *ldx, double * rcond, double *ferr, double *berr, double *work, int * iwork, int *info);
int F77_CALL(dgttrf)(int *n, double *dl, double *d__,  double *du, double *du2, int *ipiv, int *info);
int F77_CALL(dgttrs)(char *trans, int *n, int *nrhs,  double *dl, double *d__, double *du, double *du2,  int *ipiv, double *b, int *ldb, int *info);
int F77_CALL(dgtts2)(int *itrans, int *n, int *nrhs,  double *dl, double *d__, double *du, double *du2,  int *ipiv, double *b, int *ldb);
int F77_CALL(dhgeqz)(char *job, char *compq, char *compz, int *n,  int *ilo, int *ihi, double *a, int *lda, double * b, int *ldb, double *alphar, double *alphai, double * beta, double *q, int *ldq, double *z__, int *ldz,  double *work, int *lwork, int *info);
int F77_CALL(dhsein)(char *side, char *eigsrc, char *initv, logical * select, int *n, double *h__, int *ldh, double *wr,  double *wi, double *vl, int *ldvl, double *vr,  int *ldvr, int *mm, int *m, double *work, int * ifaill, int *ifailr, int *info);
int F77_CALL(dhseqr)(char *job, char *compz, int *n, int *ilo, int *ihi, double *h__, int *ldh, double *wr,  double *wi, double *z__, int *ldz, double *work,  int *lwork, int *info);
int F77_CALL(dlabad)(double *small, double *large);
int F77_CALL(dlabrd)(int *m, int *n, int *nb, double * a, int *lda, double *d__, double *e, double *tauq,  double *taup, double *x, int *ldx, double *y, int  *ldy);
int F77_CALL(dlacon)(int *n, double *v, double *x,  int *isgn, double *est, int *kase);
int F77_CALL(dlacpy)(char *uplo, int *m, int *n, double * a, int *lda, double *b, int *ldb);
int F77_CALL(dladiv)(double *a, double *b, double *c__,  double *d__, double *p, double *q);
int F77_CALL(dlae2)(double *a, double *b, double *c__,  double *rt1, double *rt2);
int F77_CALL(dlaebz)(int *ijob, int *nitmax, int *n,  int *mmax, int *minp, int *nbmin, double *abstol,  double *reltol, double *pivmin, double *d__, double * e, double *e2, int *nval, double *ab, double *c__,  int *mout, int *nab, double *work, int *iwork,  int *info);
int F77_CALL(dlaed0)(int *icompq, int *qsiz, int *n,  double *d__, double *e, double *q, int *ldq,  double *qstore, int *ldqs, double *work, int *iwork,  int *info);
int F77_CALL(dlaed1)(int *n, double *d__, double *q,  int *ldq, int *indxq, double *rho, int *cutpnt,  double *work, int *iwork, int *info);
int F77_CALL(dlaed2)(int *k, int *n, int *n1, double * d__, double *q, int *ldq, int *indxq, double *rho,  double *z__, double *dlamda, double *w, double *q2,  int *indx, int *indxc, int *indxp, int *coltyp,  int *info);
int F77_CALL(dlaed3)(int *k, int *n, int *n1, double * d__, double *q, int *ldq, double *rho, double *dlamda, double *q2, int *indx, int *ctot, double *w,  double *s, int *info);
int F77_CALL(dlaed4)(int *n, int *i__, double *d__,  double *z__, double *delta, double *rho, double *dlam, int *info);
int F77_CALL(dlaed5)(int *i__, double *d__, double *z__,  double *delta, double *rho, double *dlam);
int F77_CALL(dlaed6)(int *kniter, logical *orgati, double * rho, double *d__, double *z__, double *finit, double * tau, int *info);
int F77_CALL(dlaed7)(int *icompq, int *n, int *qsiz,  int *tlvls, int *curlvl, int *curpbm, double *d__,  double *q, int *ldq, int *indxq, double *rho, int  *cutpnt, double *qstore, int *qptr, int *prmptr, int * perm, int *givptr, int *givcol, double *givnum,  double *work, int *iwork, int *info);
int F77_CALL(dlaed8)(int *icompq, int *k, int *n, int  *qsiz, double *d__, double *q, int *ldq, int *indxq,  double *rho, int *cutpnt, double *z__, double *dlamda, double *q2, int *ldq2, double *w, int *perm, int  *givptr, int *givcol, double *givnum, int *indxp, int  *indx, int *info);
int F77_CALL(dlaed9)(int *k, int *kstart, int *kstop,  int *n, double *d__, double *q, int *ldq, double * rho, double *dlamda, double *w, double *s, int *lds,  int *info);
int F77_CALL(dlaeda)(int *n, int *tlvls, int *curlvl,  int *curpbm, int *prmptr, int *perm, int *givptr,  int *givcol, double *givnum, double *q, int *qptr,  double *z__, double *ztemp, int *info);
int F77_CALL(dlaein)(logical *rightv, logical *noinit, int *n,  double *h__, int *ldh, double *wr, double *wi,  double *vr, double *vi, double *b, int *ldb,  double *work, double *eps3, double *smlnum, double * bignum, int *info);
int F77_CALL(dlaev2)(double *a, double *b, double *c__,  double *rt1, double *rt2, double *cs1, double *sn1);
int F77_CALL(dlaexc)(logical *wantq, int *n, double *t,  int *ldt, double *q, int *ldq, int *j1, int *n1,  int *n2, double *work, int *info);
int F77_CALL(dlag2)(double *a, int *lda, double *b,  int *ldb, double *safmin, double *scale1, double * scale2, double *wr1, double *wr2, double *wi);
int F77_CALL(dlags2)(logical *upper, double *a1, double *a2,  double *a3, double *b1, double *b2, double *b3,  double *csu, double *snu, double *csv, double *snv,  double *csq, double *snq);
int F77_CALL(dlagtf)(int *n, double *a, double *lambda,  double *b, double *c__, double *tol, double *d__,  int *in, int *info);
int F77_CALL(dlagtm)(char *trans, int *n, int *nrhs,  double *alpha, double *dl, double *d__, double *du,  double *x, int *ldx, double *beta, double *b, int  *ldb);
int F77_CALL(dlagts)(int *job, int *n, double *a,  double *b, double *c__, double *d__, int *in,  double *y, double *tol, int *info);
int F77_CALL(dlagv2)(double *a, int *lda, double *b,  int *ldb, double *alphar, double *alphai, double * beta, double *csl, double *snl, double *csr, double * snr);
int F77_CALL(dlahqr)(logical *wantt, logical *wantz, int *n,  int *ilo, int *ihi, double *h__, int *ldh, double  *wr, double *wi, int *iloz, int *ihiz, double *z__,  int *ldz, int *info);
int F77_CALL(dlahrd)(int *n, int *k, int *nb, double * a, int *lda, double *tau, double *t, int *ldt,  double *y, int *ldy);
int F77_CALL(dlaic1)(int *job, int *j, double *x,  double *sest, double *w, double *gamma, double * sestpr, double *s, double *c__);
int F77_CALL(dlaln2)(logical *ltrans, int *na, int *nw,  double *smin, double *ca, double *a, int *lda,  double *d1, double *d2, double *b, int *ldb,  double *wr, double *wi, double *x, int *ldx,  double *scale, double *xnorm, int *info);
int F77_CALL(dlals0)(int *icompq, int *nl, int *nr,  int *sqre, int *nrhs, double *b, int *ldb, double  *bx, int *ldbx, int *perm, int *givptr, int *givcol,  int *ldgcol, double *givnum, int *ldgnum, double * poles, double *difl, double *difr, double *z__, int * k, double *c__, double *s, double *work, int *info);
int F77_CALL(dlalsa)(int *icompq, int *smlsiz, int *n,  int *nrhs, double *b, int *ldb, double *bx, int * ldbx, double *u, int *ldu, double *vt, int *k,  double *difl, double *difr, double *z__, double * poles, int *givptr, int *givcol, int *ldgcol, int * perm, double *givnum, double *c__, double *s, double * work, int *iwork, int *info);
int F77_CALL(dlalsd)(char *uplo, int *smlsiz, int *n, int  *nrhs, double *d__, double *e, double *b, int *ldb,  double *rcond, int *rank, double *work, int *iwork,  int *info);
int F77_CALL(dlamc1)(int *beta, int *t, logical *rnd, logical  *ieee1);
int F77_CALL(dlamc2)(int *beta, int *t, logical *rnd,  double *eps, int *emin, double *rmin, int *emax,  double *rmax);
int F77_CALL(dlamc4)(int *emin, double *start, int *base);
int F77_CALL(dlamc5)(int *beta, int *p, int *emin,  logical *ieee, int *emax, double *rmax);
int F77_CALL(dlamrg)(int *n1, int *n2, double *a, int  *dtrd1, int *dtrd2, int *index);
int F77_CALL(dlanv2)(double *a, double *b, double *c__,  double *d__, double *rt1r, double *rt1i, double *rt2r, double *rt2i, double *cs, double *sn);
int F77_CALL(dlapll)(int *n, double *x, int *incx,  double *y, int *incy, double *ssmin);
int F77_CALL(dlapmt)(logical *forwrd, int *m, int *n,  double *x, int *ldx, int *k);
int F77_CALL(dlaqgb)(int *m, int *n, int *kl, int *ku, double *ab, int *ldab, double *r__, double *c__,  double *rowcnd, double *colcnd, double *amax, char *equed);
int F77_CALL(dlaqge)(int *m, int *n, double *a, int * lda, double *r__, double *c__, double *rowcnd, double  *colcnd, double *amax, char *equed);
int F77_CALL(dlaqp2)(int *m, int *n, int *offset,  double *a, int *lda, int *jpvt, double *tau,  double *vn1, double *vn2, double *work);
int F77_CALL(dlaqps)(int *m, int *n, int *offset, int  *nb, int *kb, double *a, int *lda, int *jpvt,  double *tau, double *vn1, double *vn2, double *auxv,  double *f, int *ldf);
int F77_CALL(dlaqsb)(char *uplo, int *n, int *kd, double * ab, int *ldab, double *s, double *scond, double *amax, char *equed);
int F77_CALL(dlaqsp)(char *uplo, int *n, double *ap,  double *s, double *scond, double *amax, char *equed);
int F77_CALL(dlaqsy)(char *uplo, int *n, double *a, int * lda, double *s, double *scond, double *amax, char *equed);
int F77_CALL(dlaqtr)(logical *ltran, logical *lreal, int *n,  double *t, int *ldt, double *b, double *w, double  *scale, double *x, double *work, int *info);
int F77_CALL(dlar1v)(int *n, int *b1, int *bn, double  *sigma, double *d__, double *l, double *ld, double * lld, double *gersch, double *z__, double *ztz, double  *mingma, int *r__, int *isuppz, double *work);
int F77_CALL(dlar2v)(int *n, double *x, double *y,  double *z__, int *incx, double *c__, double *s,  int *incc);
int F77_CALL(dlarf)(char *side, int *m, int *n, double *v, int *incv, double *tau, double *c__, int *ldc,  double *work);
int F77_CALL(dlarfb)(char *side, char *trans, char *direct, char * storev, int *m, int *n, int *k, double *v, int * ldv, double *t, int *ldt, double *c__, int *ldc,  double *work, int *ldwork);
int F77_CALL(dlarfg)(int *n, double *alpha, double *x,  int *incx, double *tau);
int F77_CALL(dlarft)(char *direct, char *storev, int *n, int * k, double *v, int *ldv, double *tau, double *t,  int *ldt);
int F77_CALL(dlarfx)(char *side, int *m, int *n, double * v, double *tau, double *c__, int *ldc, double *work);
int F77_CALL(dlargv)(int *n, double *x, int *incx,  double *y, int *incy, double *c__, int *incc);
int F77_CALL(dlarnv)(int *idist, int *iseed, int *n,  double *x);
int F77_CALL(dlarrb)(int *n, double *d__, double *l,  double *ld, double *lld, int *ifirst, int *ilast,  double *sigma, double *reltol, double *w, double * wgap, double *werr, double *work, int *iwork, int * info);
int F77_CALL(dlarre)(int *n, double *d__, double *e,  double *tol, int *nsplit, int *isplit, int *m,  double *w, double *woff, double *gersch, double *work, int *info);
int F77_CALL(dlarrf)(int *n, double *d__, double *l,  double *ld, double *lld, int *ifirst, int *ilast,  double *w, double *dplus, double *lplus, double *work, int *iwork, int *info);
int F77_CALL(dlarrv)(int *n, double *d__, double *l,  int *isplit, int *m, double *w, int *iblock,  double *gersch, double *tol, double *z__, int *ldz,  int *isuppz, double *work, int *iwork, int *info);
int F77_CALL(dlartg)(double *f, double *g, double *cs,  double *sn, double *r__);
int F77_CALL(dlartv)(int *n, double *x, int *incx,  double *y, int *incy, double *c__, double *s, int  *incc);
int F77_CALL(dlaruv)(int *iseed, int *n, double *x);
int F77_CALL(dlarz)(char *side, int *m, int *n, int *l,  double *v, int *incv, double *tau, double *c__,  int *ldc, double *work);
int F77_CALL(dlarzb)(char *side, char *trans, char *direct, char * storev, int *m, int *n, int *k, int *l, double *v, int *ldv, double *t, int *ldt, double *c__, int * ldc, double *work, int *ldwork);
int F77_CALL(dlarzt)(char *direct, char *storev, int *n, int * k, double *v, int *ldv, double *tau, double *t,  int *ldt);
int F77_CALL(dlas2)(double *f, double *g, double *h__,  double *ssmin, double *ssmax);
int F77_CALL(dlascl)(char *type__, int *kl, int *ku,  double *cfrom, double *cto, int *m, int *n,  double *a, int *lda, int *info);
int F77_CALL(dlasd0)(int *n, int *sqre, double *d__,  double *e, double *u, int *ldu, double *vt, int * ldvt, int *smlsiz, int *iwork, double *work, int * info);
int F77_CALL(dlasd1)(int *nl, int *nr, int *sqre,  double *d__, double *alpha, double *beta, double *u,  int *ldu, double *vt, int *ldvt, int *idxq, int * iwork, double *work, int *info);
int F77_CALL(dlasd2)(int *nl, int *nr, int *sqre, int  *k, double *d__, double *z__, double *alpha, double * beta, double *u, int *ldu, double *vt, int *ldvt,  double *dsigma, double *u2, int *ldu2, double *vt2,  int *ldvt2, int *idxp, int *idx, int *idxc, int * idxq, int *coltyp, int *info);
int F77_CALL(dlasd3)(int *nl, int *nr, int *sqre, int  *k, double *d__, double *q, int *ldq, double *dsigma,  double *u, int *ldu, double *u2, int *ldu2,  double *vt, int *ldvt, double *vt2, int *ldvt2,  int *idxc, int *ctot, double *z__, int *info);
int F77_CALL(dlasd4)(int *n, int *i__, double *d__,  double *z__, double *delta, double *rho, double * sigma, double *work, int *info);
int F77_CALL(dlasd5)(int *i__, double *d__, double *z__,  double *delta, double *rho, double *dsigma, double * work);
int F77_CALL(dlasd6)(int *icompq, int *nl, int *nr,  int *sqre, double *d__, double *vf, double *vl,  double *alpha, double *beta, int *idxq, int *perm,  int *givptr, int *givcol, int *ldgcol, double *givnum, int *ldgnum, double *poles, double *difl, double * difr, double *z__, int *k, double *c__, double *s,  double *work, int *iwork, int *info);
int F77_CALL(dlasd7)(int *icompq, int *nl, int *nr,  int *sqre, int *k, double *d__, double *z__,  double *zw, double *vf, double *vfw, double *vl,  double *vlw, double *alpha, double *beta, double * dsigma, int *idx, int *idxp, int *idxq, int *perm,  int *givptr, int *givcol, int *ldgcol, double *givnum, int *ldgnum, double *c__, double *s, int *info);
int F77_CALL(dlasd8)(int *icompq, int *k, double *d__,  double *z__, double *vf, double *vl, double *difl,  double *difr, int *lddifr, double *dsigma, double * work, int *info);
int F77_CALL(dlasd9)(int *icompq, int *ldu, int *k,  double *d__, double *z__, double *vf, double *vl,  double *difl, double *difr, double *dsigma, double * work, int *info);
int F77_CALL(dlasda)(int *icompq, int *smlsiz, int *n,  int *sqre, double *d__, double *e, double *u, int  *ldu, double *vt, int *k, double *difl, double *difr,  double *z__, double *poles, int *givptr, int *givcol,  int *ldgcol, int *perm, double *givnum, double *c__,  double *s, double *work, int *iwork, int *info);
int F77_CALL(dlasdq)(char *uplo, int *sqre, int *n, int * ncvt, int *nru, int *ncc, double *d__, double *e,  double *vt, int *ldvt, double *u, int *ldu,  double *c__, int *ldc, double *work, int *info);
int F77_CALL(dlasdt)(int *n, int *lvl, int *nd, int * inode, int *ndiml, int *ndimr, int *msub);
int F77_CALL(dlaset)(char *uplo, int *m, int *n, double * alpha, double *beta, double *a, int *lda);
int F77_CALL(dlasq1)(int *n, double *d__, double *e,  double *work, int *info);
int F77_CALL(dlasq2)(int *n, double *z__, int *info);
int F77_CALL(dlasq3)(int *i0, int *n0, double *z__,  int *pp, double *dmin__, double *sigma, double *desig, double *qmax, int *nfail, int *iter, int *ndiv,  logical *ieee);
int F77_CALL(dlasq4)(int *i0, int *n0, double *z__,  int *pp, int *n0in, double *dmin__, double *dmin1,  double *dmin2, double *dn, double *dn1, double *dn2,  double *tau, int *ttype);
int F77_CALL(dlasq5)(int *i0, int *n0, double *z__,  int *pp, double *tau, double *dmin__, double *dmin1,  double *dmin2, double *dn, double *dnm1, double *dnm2, logical *ieee);
int F77_CALL(dlasq6)(int *i0, int *n0, double *z__,  int *pp, double *dmin__, double *dmin1, double *dmin2, double *dn, double *dnm1, double *dnm2);
int F77_CALL(dlasr)(char *side, char *pivot, char *direct, int *m, int *n, double *c__, double *s, double *a, int * lda);
int F77_CALL(dlasrt)(char *id, int *n, double *d__, int * info);
int F77_CALL(dlassq)(int *n, double *x, int *incx,  double *scale, double *sumsq);
int F77_CALL(dlasv2)(double *f, double *g, double *h__,  double *ssmin, double *ssmax, double *snr, double * csr, double *snl, double *csl);
int F77_CALL(dlaswp)(int *n, double *a, int *lda, int  *k1, int *k2, int *ipiv, int *incx);
int F77_CALL(dlasy2)(logical *ltranl, logical *ltranr, int *isgn,  int *n1, int *n2, double *tl, int *ldtl, double * tr, int *ldtr, double *b, int *ldb, double *scale,  double *x, int *ldx, double *xnorm, int *info);
int F77_CALL(dlasyf)(char *uplo, int *n, int *nb, int *kb, double *a, int *lda, int *ipiv, double *w, int * ldw, int *info);
int F77_CALL(dlatbs)(char *uplo, char *trans, char *diag, char * normin, int *n, int *kd, double *ab, int *ldab,  double *x, double *scale, double *cnorm, int *info);
int F77_CALL(dlatdf)(int *ijob, int *n, double *z__,  int *ldz, double *rhs, double *rdsum, double *rdscal,  int *ipiv, int *jpiv);
int F77_CALL(dlatps)(char *uplo, char *trans, char *diag, char * normin, int *n, double *ap, double *x, double *scale,  double *cnorm, int *info);
int F77_CALL(dlatrd)(char *uplo, int *n, int *nb, double * a, int *lda, double *e, double *tau, double *w,  int *ldw);
int F77_CALL(dlatrs)(char *uplo, char *trans, char *diag, char * normin, int *n, double *a, int *lda, double *x,  double *scale, double *cnorm, int *info);
int F77_CALL(dlatrz)(int *m, int *n, int *l, double * a, int *lda, double *tau, double *work);
int F77_CALL(dlatzm)(char *side, int *m, int *n, double * v, int *incv, double *tau, double *c1, double *c2,  int *ldc, double *work);
int F77_CALL(dlauu2)(char *uplo, int *n, double *a, int * lda, int *info);
int F77_CALL(dlauum)(char *uplo, int *n, double *a, int * lda, int *info);
int F77_CALL(dopgtr)(char *uplo, int *n, double *ap,  double *tau, double *q, int *ldq, double *work,  int *info);
int F77_CALL(dopmtr)(char *side, char *uplo, char *trans, int *m,  int *n, double *ap, double *tau, double *c__, int  *ldc, double *work, int *info);
int F77_CALL(dorg2l)(int *m, int *n, int *k, double * a, int *lda, double *tau, double *work, int *info);
int F77_CALL(dorg2r)(int *m, int *n, int *k, double * a, int *lda, double *tau, double *work, int *info);
int F77_CALL(dorgbr)(char *vect, int *m, int *n, int *k,  double *a, int *lda, double *tau, double *work,  int *lwork, int *info);
int F77_CALL(dorghr)(int *n, int *ilo, int *ihi,  double *a, int *lda, double *tau, double *work,  int *lwork, int *info);
int F77_CALL(dorgl2)(int *m, int *n, int *k, double * a, int *lda, double *tau, double *work, int *info);
int F77_CALL(dorglq)(int *m, int *n, int *k, double * a, int *lda, double *tau, double *work, int *lwork,  int *info);
int F77_CALL(dorgql)(int *m, int *n, int *k, double * a, int *lda, double *tau, double *work, int *lwork,  int *info);
int F77_CALL(dorgqr)(int *m, int *n, int *k, double * a, int *lda, double *tau, double *work, int *lwork,  int *info);
int F77_CALL(dorgr2)(int *m, int *n, int *k, double * a, int *lda, double *tau, double *work, int *info);
int F77_CALL(dorgrq)(int *m, int *n, int *k, double * a, int *lda, double *tau, double *work, int *lwork,  int *info);
int F77_CALL(dorgtr)(char *uplo, int *n, double *a, int * lda, double *tau, double *work, int *lwork, int *info);
int F77_CALL(dorm2l)(char *side, char *trans, int *m, int *n,  int *k, double *a, int *lda, double *tau, double * c__, int *ldc, double *work, int *info);
int F77_CALL(dorm2r)(char *side, char *trans, int *m, int *n,  int *k, double *a, int *lda, double *tau, double * c__, int *ldc, double *work, int *info);
int F77_CALL(dormbr)(char *vect, char *side, char *trans, int *m,  int *n, int *k, double *a, int *lda, double *tau,  double *c__, int *ldc, double *work, int *lwork,  int *info);
int F77_CALL(dormhr)(char *side, char *trans, int *m, int *n,  int *ilo, int *ihi, double *a, int *lda, double * tau, double *c__, int *ldc, double *work, int *lwork,  int *info);
int F77_CALL(dorml2)(char *side, char *trans, int *m, int *n,  int *k, double *a, int *lda, double *tau, double * c__, int *ldc, double *work, int *info);
int F77_CALL(dormlq)(char *side, char *trans, int *m, int *n,  int *k, double *a, int *lda, double *tau, double * c__, int *ldc, double *work, int *lwork, int *info);
int F77_CALL(dormql)(char *side, char *trans, int *m, int *n,  int *k, double *a, int *lda, double *tau, double * c__, int *ldc, double *work, int *lwork, int *info);
int F77_CALL(dormqr)(char *side, char *trans, int *m, int *n,  int *k, double *a, int *lda, double *tau, double * c__, int *ldc, double *work, int *lwork, int *info);
int F77_CALL(dormr2)(char *side, char *trans, int *m, int *n,  int *k, double *a, int *lda, double *tau, double * c__, int *ldc, double *work, int *info);
int F77_CALL(dormr3)(char *side, char *trans, int *m, int *n,  int *k, int *l, double *a, int *lda, double *tau,  double *c__, int *ldc, double *work, int *info);
int F77_CALL(dormrq)(char *side, char *trans, int *m, int *n,  int *k, double *a, int *lda, double *tau, double * c__, int *ldc, double *work, int *lwork, int *info);
int F77_CALL(dormrz)(char *side, char *trans, int *m, int *n,  int *k, int *l, double *a, int *lda, double *tau,  double *c__, int *ldc, double *work, int *lwork,  int *info);
int F77_CALL(dormtr)(char *side, char *uplo, char *trans, int *m,  int *n, double *a, int *lda, double *tau, double * c__, int *ldc, double *work, int *lwork, int *info);
int F77_CALL(dpbcon)(char *uplo, int *n, int *kd, double * ab, int *ldab, double *anorm, double *rcond, double * work, int *iwork, int *info);
int F77_CALL(dpbequ)(char *uplo, int *n, int *kd, double * ab, int *ldab, double *s, double *scond, double *amax, int *info);
int F77_CALL(dpbrfs)(char *uplo, int *n, int *kd, int * nrhs, double *ab, int *ldab, double *afb, int *ldafb,  double *b, int *ldb, double *x, int *ldx, double * ferr, double *berr, double *work, int *iwork, int * info);
int F77_CALL(dpbstf)(char *uplo, int *n, int *kd, double * ab, int *ldab, int *info);
int F77_CALL(dpbsv)(char *uplo, int *n, int *kd, int * nrhs, double *ab, int *ldab, double *b, int *ldb,  int *info);
int F77_CALL(dpbsvx)(char *fact, char *uplo, int *n, int *kd,  int *nrhs, double *ab, int *ldab, double *afb,  int *ldafb, char *equed, double *s, double *b, int * ldb, double *x, int *ldx, double *rcond, double *ferr, double *berr, double *work, int *iwork, int *info);
int F77_CALL(dpbtf2)(char *uplo, int *n, int *kd, double * ab, int *ldab, int *info);
int F77_CALL(dpbtrf)(char *uplo, int *n, int *kd, double * ab, int *ldab, int *info);
int F77_CALL(dpbtrs)(char *uplo, int *n, int *kd, int * nrhs, double *ab, int *ldab, double *b, int *ldb,  int *info);
int F77_CALL(dpocon)(char *uplo, int *n, double *a, int * lda, double *anorm, double *rcond, double *work, int * iwork, int *info);
int F77_CALL(dpoequ)(int *n, double *a, int *lda,  double *s, double *scond, double *amax, int *info);
int F77_CALL(dporfs)(char *uplo, int *n, int *nrhs,  double *a, int *lda, double *af, int *ldaf,  double *b, int *ldb, double *x, int *ldx, double * ferr, double *berr, double *work, int *iwork, int * info);
int F77_CALL(dposv)(char *uplo, int *n, int *nrhs, double  *a, int *lda, double *b, int *ldb, int *info);
int F77_CALL(dposvx)(char *fact, char *uplo, int *n, int * nrhs, double *a, int *lda, double *af, int *ldaf,  char *equed, double *s, double *b, int *ldb, double * x, int *ldx, double *rcond, double *ferr, double * berr, double *work, int *iwork, int *info);
int F77_CALL(dpotf2)(char *uplo, int *n, double *a, int * lda, int *info);
int F77_CALL(dpotrf)(char *uplo, int *n, double *a, int * lda, int *info);
int F77_CALL(dpotri)(char *uplo, int *n, double *a, int * lda, int *info);
int F77_CALL(dpotrs)(char *uplo, int *n, int *nrhs,  double *a, int *lda, double *b, int *ldb, int * info);
int F77_CALL(dppcon)(char *uplo, int *n, double *ap,  double *anorm, double *rcond, double *work, int * iwork, int *info);
int F77_CALL(dppequ)(char *uplo, int *n, double *ap,  double *s, double *scond, double *amax, int *info);
int F77_CALL(dpprfs)(char *uplo, int *n, int *nrhs,  double *ap, double *afp, double *b, int *ldb,  double *x, int *ldx, double *ferr, double *berr,  double *work, int *iwork, int *info);
int F77_CALL(dppsv)(char *uplo, int *n, int *nrhs, double  *ap, double *b, int *ldb, int *info);
int F77_CALL(dppsvx)(char *fact, char *uplo, int *n, int * nrhs, double *ap, double *afp, char *equed, double *s,  double *b, int *ldb, double *x, int *ldx, double * rcond, double *ferr, double *berr, double *work, int * iwork, int *info);
int F77_CALL(dpptrf)(char *uplo, int *n, double *ap, int * info);
int F77_CALL(dpptri)(char *uplo, int *n, double *ap, int * info);
int F77_CALL(dpptrs)(char *uplo, int *n, int *nrhs,  double *ap, double *b, int *ldb, int *info);
int F77_CALL(dptcon)(int *n, double *d__, double *e,  double *anorm, double *rcond, double *work, int *info);
int F77_CALL(dpteqr)(char *compz, int *n, double *d__,  double *e, double *z__, int *ldz, double *work,  int *info);
int F77_CALL(dptrfs)(int *n, int *nrhs, double *d__,  double *e, double *df, double *ef, double *b, int  *ldb, double *x, int *ldx, double *ferr, double *berr, double *work, int *info);
int F77_CALL(dptsv)(int *n, int *nrhs, double *d__,  double *e, double *b, int *ldb, int *info);
int F77_CALL(dptsvx)(char *fact, int *n, int *nrhs,  double *d__, double *e, double *df, double *ef,  double *b, int *ldb, double *x, int *ldx, double * rcond, double *ferr, double *berr, double *work, int * info);
int F77_CALL(dpttrf)(int *n, double *d__, double *e,  int *info);
int F77_CALL(dpttrs)(int *n, int *nrhs, double *d__,  double *e, double *b, int *ldb, int *info);
int F77_CALL(dptts2)(int *n, int *nrhs, double *d__,  double *e, double *b, int *ldb);
int F77_CALL(drscl)(int *n, double *sa, double *sx,  int *incx);
int F77_CALL(dsbev)(char *jobz, char *uplo, int *n, int *kd,  double *ab, int *ldab, double *w, double *z__,  int *ldz, double *work, int *info);
int F77_CALL(dsbevd)(char *jobz, char *uplo, int *n, int *kd,  double *ab, int *ldab, double *w, double *z__,  int *ldz, double *work, int *lwork, int *iwork,  int *liwork, int *info);
int F77_CALL(dsbevx)(char *jobz, char *range, char *uplo, int *n,  int *kd, double *ab, int *ldab, double *q, int * ldq, double *vl, double *vu, int *il, int *iu,  double *abstol, int *m, double *w, double *z__,  int *ldz, double *work, int *iwork, int *ifail,  int *info);
int F77_CALL(dsbgst)(char *vect, char *uplo, int *n, int *ka,  int *kb, double *ab, int *ldab, double *bb, int * ldbb, double *x, int *ldx, double *work, int *info);
int F77_CALL(dsbgv)(char *jobz, char *uplo, int *n, int *ka,  int *kb, double *ab, int *ldab, double *bb, int * ldbb, double *w, double *z__, int *ldz, double *work,  int *info);
int F77_CALL(dsbgvd)(char *jobz, char *uplo, int *n, int *ka,  int *kb, double *ab, int *ldab, double *bb, int * ldbb, double *w, double *z__, int *ldz, double *work,  int *lwork, int *iwork, int *liwork, int *info);
int F77_CALL(dsbgvx)(char *jobz, char *range, char *uplo, int *n,  int *ka, int *kb, double *ab, int *ldab, double * bb, int *ldbb, double *q, int *ldq, double *vl,  double *vu, int *il, int *iu, double *abstol, int  *m, double *w, double *z__, int *ldz, double *work,  int *iwork, int *ifail, int *info);
int F77_CALL(dsbtrd)(char *vect, char *uplo, int *n, int *kd,  double *ab, int *ldab, double *d__, double *e,  double *q, int *ldq, double *work, int *info);
int F77_CALL(dspcon)(char *uplo, int *n, double *ap, int * ipiv, double *anorm, double *rcond, double *work, int  *iwork, int *info);
int F77_CALL(dspev)(char *jobz, char *uplo, int *n, double * ap, double *w, double *z__, int *ldz, double *work,  int *info);
int F77_CALL(dspevd)(char *jobz, char *uplo, int *n, double * ap, double *w, double *z__, int *ldz, double *work,  int *lwork, int *iwork, int *liwork, int *info);
int F77_CALL(dspevx)(char *jobz, char *range, char *uplo, int *n,  double *ap, double *vl, double *vu, int *il, int * iu, double *abstol, int *m, double *w, double *z__,  int *ldz, double *work, int *iwork, int *ifail,  int *info);
int F77_CALL(dspgst)(int *itype, char *uplo, int *n,  double *ap, double *bp, int *info);
int F77_CALL(dspgv)(int *itype, char *jobz, char *uplo, int * n, double *ap, double *bp, double *w, double *z__,  int *ldz, double *work, int *info);
int F77_CALL(dspgvd)(int *itype, char *jobz, char *uplo, int * n, double *ap, double *bp, double *w, double *z__,  int *ldz, double *work, int *lwork, int *iwork,  int *liwork, int *info);
int F77_CALL(dspgvx)(int *itype, char *jobz, char *range, char * uplo, int *n, double *ap, double *bp, double *vl,  double *vu, int *il, int *iu, double *abstol, int  *m, double *w, double *z__, int *ldz, double *work,  int *iwork, int *ifail, int *info);
int F77_CALL(dsprfs)(char *uplo, int *n, int *nrhs,  double *ap, double *afp, int *ipiv, double *b,  int *ldb, double *x, int *ldx, double *ferr,  double *berr, double *work, int *iwork, int *info);
int F77_CALL(dspsv)(char *uplo, int *n, int *nrhs, double  *ap, int *ipiv, double *b, int *ldb, int *info);
int F77_CALL(dspsvx)(char *fact, char *uplo, int *n, int * nrhs, double *ap, double *afp, int *ipiv, double *b,  int *ldb, double *x, int *ldx, double *rcond,  double *ferr, double *berr, double *work, int *iwork,  int *info);
int F77_CALL(dsptrd)(char *uplo, int *n, double *ap,  double *d__, double *e, double *tau, int *info);
int F77_CALL(dsptrf)(char *uplo, int *n, double *ap, int * ipiv, int *info);
int F77_CALL(dsptri)(char *uplo, int *n, double *ap, int * ipiv, double *work, int *info);
int F77_CALL(dsptrs)(char *uplo, int *n, int *nrhs,  double *ap, int *ipiv, double *b, int *ldb, int * info);
int F77_CALL(dstebz)(char *range, char *order, int *n, double  *vl, double *vu, int *il, int *iu, double *abstol,  double *d__, double *e, int *m, int *nsplit,  double *w, int *iblock, int *isplit, double *work,  int *iwork, int *info);
int F77_CALL(dstedc)(char *compz, int *n, double *d__,  double *e, double *z__, int *ldz, double *work,  int *lwork, int *iwork, int *liwork, int *info);
int F77_CALL(dstegr)(char *jobz, char *range, int *n, double * d__, double *e, double *vl, double *vu, int *il,  int *iu, double *abstol, int *m, double *w,  double *z__, int *ldz, int *isuppz, double *work,  int *lwork, int *iwork, int *liwork, int *info);
int F77_CALL(dstein)(int *n, double *d__, double *e,  int *m, double *w, int *iblock, int *isplit,  double *z__, int *ldz, double *work, int *iwork,  int *ifail, int *info);
int F77_CALL(dsteqr)(char *compz, int *n, double *d__,  double *e, double *z__, int *ldz, double *work,  int *info);
int F77_CALL(dsterf)(int *n, double *d__, double *e,  int *info);
int F77_CALL(dstev)(char *jobz, int *n, double *d__,  double *e, double *z__, int *ldz, double *work,  int *info);
int F77_CALL(dstevd)(char *jobz, int *n, double *d__,  double *e, double *z__, int *ldz, double *work,  int *lwork, int *iwork, int *liwork, int *info);
int F77_CALL(dstevr)(char *jobz, char *range, int *n, double * d__, double *e, double *vl, double *vu, int *il,  int *iu, double *abstol, int *m, double *w,  double *z__, int *ldz, int *isuppz, double *work,  int *lwork, int *iwork, int *liwork, int *info);
int F77_CALL(dstevx)(char *jobz, char *range, int *n, double * d__, double *e, double *vl, double *vu, int *il,  int *iu, double *abstol, int *m, double *w,  double *z__, int *ldz, double *work, int *iwork,  int *ifail, int *info);
int F77_CALL(dsycon)(char *uplo, int *n, double *a, int * lda, int *ipiv, double *anorm, double *rcond, double * work, int *iwork, int *info);
int F77_CALL(dsyev)(char *jobz, char *uplo, int *n, double *a, int *lda, double *w, double *work, int *lwork,  int *info);
int F77_CALL(dsyevd)(char *jobz, char *uplo, int *n, double * a, int *lda, double *w, double *work, int *lwork,  int *iwork, int *liwork, int *info);
int F77_CALL(dsyevr)(char *jobz, char *range, char *uplo, int *n,  double *a, int *lda, double *vl, double *vu, int * il, int *iu, double *abstol, int *m, double *w,  double *z__, int *ldz, int *isuppz, double *work,  int *lwork, int *iwork, int *liwork, int *info);
int F77_CALL(dsyevx)(char *jobz, char *range, char *uplo, int *n,  double *a, int *lda, double *vl, double *vu, int * il, int *iu, double *abstol, int *m, double *w,  double *z__, int *ldz, double *work, int *lwork,  int *iwork, int *ifail, int *info);
int F77_CALL(dsygs2)(int *itype, char *uplo, int *n,  double *a, int *lda, double *b, int *ldb, int * info);
int F77_CALL(dsygst)(int *itype, char *uplo, int *n,  double *a, int *lda, double *b, int *ldb, int * info);
int F77_CALL(dsygv)(int *itype, char *jobz, char *uplo, int * n, double *a, int *lda, double *b, int *ldb,  double *w, double *work, int *lwork, int *info);
int F77_CALL(dsygvd)(int *itype, char *jobz, char *uplo, int * n, double *a, int *lda, double *b, int *ldb,  double *w, double *work, int *lwork, int *iwork,  int *liwork, int *info);
int F77_CALL(dsygvx)(int *itype, char *jobz, char *range, char * uplo, int *n, double *a, int *lda, double *b, int  *ldb, double *vl, double *vu, int *il, int *iu,  double *abstol, int *m, double *w, double *z__,  int *ldz, double *work, int *lwork, int *iwork,  int *ifail, int *info);
int F77_CALL(dsyrfs)(char *uplo, int *n, int *nrhs,  double *a, int *lda, double *af, int *ldaf, int * ipiv, double *b, int *ldb, double *x, int *ldx,  double *ferr, double *berr, double *work, int *iwork,  int *info);
int F77_CALL(dsysv)(char *uplo, int *n, int *nrhs, double  *a, int *lda, int *ipiv, double *b, int *ldb,  double *work, int *lwork, int *info);
int F77_CALL(dsysvx)(char *fact, char *uplo, int *n, int * nrhs, double *a, int *lda, double *af, int *ldaf,  int *ipiv, double *b, int *ldb, double *x, int * ldx, double *rcond, double *ferr, double *berr,  double *work, int *lwork, int *iwork, int *info);
int F77_CALL(dsytd2)(char *uplo, int *n, double *a, int * lda, double *d__, double *e, double *tau, int *info);
int F77_CALL(dsytf2)(char *uplo, int *n, double *a, int * lda, int *ipiv, int *info);
int F77_CALL(dsytrd)(char *uplo, int *n, double *a, int * lda, double *d__, double *e, double *tau, double * work, int *lwork, int *info);
int F77_CALL(dsytrf)(char *uplo, int *n, double *a, int * lda, int *ipiv, double *work, int *lwork, int *info);
int F77_CALL(dsytri)(char *uplo, int *n, double *a, int * lda, int *ipiv, double *work, int *info);
int F77_CALL(dsytrs)(char *uplo, int *n, int *nrhs,  double *a, int *lda, int *ipiv, double *b, int * ldb, int *info);
int F77_CALL(dtbcon)(char *norm, char *uplo, char *diag, int *n,  int *kd, double *ab, int *ldab, double *rcond,  double *work, int *iwork, int *info);
int F77_CALL(dtbrfs)(char *uplo, char *trans, char *diag, int *n,  int *kd, int *nrhs, double *ab, int *ldab, double  *b, int *ldb, double *x, int *ldx, double *ferr,  double *berr, double *work, int *iwork, int *info);
int F77_CALL(dtbtrs)(char *uplo, char *trans, char *diag, int *n,  int *kd, int *nrhs, double *ab, int *ldab, double  *b, int *ldb, int *info);
int F77_CALL(dtgevc)(char *side, char *howmny, logical *select,  int *n, double *a, int *lda, double *b, int *ldb,  double *vl, int *ldvl, double *vr, int *ldvr, int  *mm, int *m, double *work, int *info);
int F77_CALL(dtgex2)(logical *wantq, logical *wantz, int *n,  double *a, int *lda, double *b, int *ldb, double * q, int *ldq, double *z__, int *ldz, int *j1, int * n1, int *n2, double *work, int *lwork, int *info);
int F77_CALL(dtgexc)(logical *wantq, logical *wantz, int *n,  double *a, int *lda, double *b, int *ldb, double * q, int *ldq, double *z__, int *ldz, int *ifst,  int *ilst, double *work, int *lwork, int *info);
int F77_CALL(dtgsen)(int *ijob, logical *wantq, logical *wantz,  logical *select, int *n, double *a, int *lda, double * b, int *ldb, double *alphar, double *alphai, double * beta, double *q, int *ldq, double *z__, int *ldz,  int *m, double *pl, double *pr, double *dif,  double *work, int *lwork, int *iwork, int *liwork,  int *info);
int F77_CALL(dtgsja)(char *jobu, char *jobv, char *jobq, int *m,  int *p, int *n, int *k, int *l, double *a,  int *lda, double *b, int *ldb, double *tola,  double *tolb, double *alpha, double *beta, double *u,  int *ldu, double *v, int *ldv, double *q, int * ldq, double *work, int *ncycle, int *info);
int F77_CALL(dtgsna)(char *job, char *howmny, logical *select,  int *n, double *a, int *lda, double *b, int *ldb,  double *vl, int *ldvl, double *vr, int *ldvr,  double *s, double *dif, int *mm, int *m, double * work, int *lwork, int *iwork, int *info);
int F77_CALL(dtgsy2)(char *trans, int *ijob, int *m, int * n, double *a, int *lda, double *b, int *ldb,  double *c__, int *ldc, double *d__, int *ldd,  double *e, int *lde, double *f, int *ldf, double * scale, double *rdsum, double *rdscal, int *iwork, int  *pq, int *info);
int F77_CALL(dtgsyl)(char *trans, int *ijob, int *m, int * n, double *a, int *lda, double *b, int *ldb,  double *c__, int *ldc, double *d__, int *ldd,  double *e, int *lde, double *f, int *ldf, double * scale, double *dif, double *work, int *lwork, int * iwork, int *info);
int F77_CALL(dtpcon)(char *norm, char *uplo, char *diag, int *n,  double *ap, double *rcond, double *work, int *iwork,  int *info);
int F77_CALL(dtprfs)(char *uplo, char *trans, char *diag, int *n,  int *nrhs, double *ap, double *b, int *ldb,  double *x, int *ldx, double *ferr, double *berr,  double *work, int *iwork, int *info);
int F77_CALL(dtptri)(char *uplo, char *diag, int *n, double * ap, int *info);
int F77_CALL(dtptrs)(char *uplo, char *trans, char *diag, int *n,  int *nrhs, double *ap, double *b, int *ldb, int * info);
int F77_CALL(dtrcon)(char *norm, char *uplo, char *diag, int *n,  double *a, int *lda, double *rcond, double *work,  int *iwork, int *info);
int F77_CALL(dtrevc)(char *side, char *howmny, logical *select,  int *n, double *t, int *ldt, double *vl, int * ldvl, double *vr, int *ldvr, int *mm, int *m,  double *work, int *info);
int F77_CALL(dtrexc)(char *compq, int *n, double *t, int * ldt, double *q, int *ldq, int *ifst, int *ilst,  double *work, int *info);
int F77_CALL(dtrrfs)(char *uplo, char *trans, char *diag, int *n,  int *nrhs, double *a, int *lda, double *b, int * ldb, double *x, int *ldx, double *ferr, double *berr,  double *work, int *iwork, int *info);
int F77_CALL(dtrsen)(char *job, char *compq, logical *select, int  *n, double *t, int *ldt, double *q, int *ldq,  double *wr, double *wi, int *m, double *s, double  *sep, double *work, int *lwork, int *iwork, int * liwork, int *info);
int F77_CALL(dtrsna)(char *job, char *howmny, logical *select,  int *n, double *t, int *ldt, double *vl, int * ldvl, double *vr, int *ldvr, double *s, double *sep,  int *mm, int *m, double *work, int *ldwork, int * iwork, int *info);
int F77_CALL(dtrsyl)(char *trana, char *tranb, int *isgn, int  *m, int *n, double *a, int *lda, double *b, int * ldb, double *c__, int *ldc, double *scale, int *info);
int F77_CALL(dtrti2)(char *uplo, char *diag, int *n, double * a, int *lda, int *info);
int F77_CALL(dtrtri)(char *uplo, char *diag, int *n, double * a, int *lda, int *info);
int F77_CALL(dtrtrs)(char *uplo, char *trans, char *diag, int *n,  int *nrhs, double *a, int *lda, double *b, int * ldb, int *info);
int F77_CALL(dtzrqf)(int *m, int *n, double *a, int * lda, double *tau, int *info);
int F77_CALL(dtzrzf)(int *m, int *n, double *a, int * lda, double *tau, double *work, int *lwork, int *info);
int F77_CALL(ilaenv)(int *ispec, char *name__, char *opts, int *n1,  int *n2, int *n3, int *n4, ftnlen name_len, ftnlen  opts_len);
int F77_CALL(slamc1)(int *beta, int *t, logical *rnd, logical  *ieee1);
int F77_CALL(slasdt)(int *n, int *lvl, int *nd, int * inode, int *ndiml, int *ndimr, int *msub);
int F77_CALL(xerbla)(char *srname, int *info);

}


#endif /* __CLAPACK_H */
