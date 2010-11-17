
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <math.h>
#include <float.h>
#include <glib.h>
#include <check.h>
#include <checkutils.h>
#include <libla/ieee754.h>

START_TEST (test_nextup_nan)
{
    ck_assert_feq(la_nextup(la_mknan(1234)), la_mknan(1234));
}
END_TEST

START_TEST (test_nextup_neg)
{
    ck_assert_feq(la_nextup(-INFINITY), -DBL_MAX);
    ck_assert_feq(la_nextup(-1 - DBL_EPSILON), -1);
    ck_assert_feq(la_nextup(-2), -2 + DBL_EPSILON);
}
END_TEST

START_TEST (test_nextup_neg_denorm)
{
    ck_assert_feq(la_nextup(-DBL_MIN), -DBL_MIN*(1 - DBL_EPSILON));
    ck_assert_feq(la_nextup(-DBL_MIN*(1-DBL_EPSILON)),
                            -DBL_MIN*(1-2*DBL_EPSILON));
    ck_assert_feq(la_nextup(-DBL_MIN*DBL_EPSILON), -0.0);
}
END_TEST

START_TEST (test_nextup_zero)
{
    ck_assert_feq(la_nextup(-0.0), DBL_MIN * DBL_EPSILON);
    ck_assert_feq(la_nextup(+0.0), DBL_MIN * DBL_EPSILON);
}
END_TEST

START_TEST (test_nextup_pos_denorm)
{
    ck_assert_feq(la_nextup(DBL_MIN*(1-DBL_EPSILON)), DBL_MIN);
    ck_assert_feq(la_nextup(DBL_MIN), DBL_MIN*(1+DBL_EPSILON));
}
END_TEST

START_TEST (test_nextup_pos)
{
    ck_assert_feq(la_nextup(1),1 + DBL_EPSILON);
    ck_assert_feq(la_nextup(2 - DBL_EPSILON), 2);
    ck_assert_feq(la_nextup(DBL_MAX), INFINITY);
}
END_TEST


TCase *
nextup_tcase ()
{
    TCase *tc = tcase_create("nextup");
    tcase_add_test(tc, test_nextup_nan);
    tcase_add_test(tc, test_nextup_neg);
    tcase_add_test(tc, test_nextup_neg_denorm);
    tcase_add_test(tc, test_nextup_zero);
    tcase_add_test(tc, test_nextup_pos_denorm);
    tcase_add_test(tc, test_nextup_pos);
    return tc;
}

START_TEST (test_nextdown)
{
    ck_assert_feq(la_nextdown(1 + DBL_EPSILON), 1);
}
END_TEST

TCase *
nextdown_tcase ()
{
    TCase *tc = tcase_create("nextdown");
    tcase_add_test(tc, test_nextdown);
    return tc;
}


START_TEST (test_ieeemean)
{
    ck_assert(la_ieeemean(-0.0, -1e-20) < 0);
    ck_assert(la_ieeemean(+0.0, +1e-20) > 0);
    ck_assert_feq(la_ieeemean(1, 4), 2);
    ck_assert_feq(la_ieeemean(2*1.013, 8*1.013), 4*1.013);
    ck_assert_feq(la_ieeemean(-1, -4), -2);
    ck_assert_feq(la_ieeemean(-1, -2), -1.5);
    ck_assert_feq(la_ieeemean(-1*(1+8*DBL_EPSILON), -2*(1+8*DBL_EPSILON)),
                  -1.5*(1+5*DBL_EPSILON));
    ck_assert_feq(la_ieeemean(0, INFINITY), 1.5);
}
END_TEST

TCase *
ieeemean_tcase ()
{
    TCase *tc = tcase_create("ieeemean");
    tcase_add_test(tc, test_ieeemean);
    return tc;
}


/* feqrel tests ported from Tango v 0.99.9,
 * public domain code written by Don Clugston
 */

/* exact equality */
START_TEST (test_feqrel_exact)
{
	ck_assert_int_eq(la_feqrel(DBL_MAX, DBL_MAX), DBL_MANT_DIG);
	ck_assert_int_eq(la_feqrel(0.0, 0.0), DBL_MANT_DIG);
	ck_assert_int_eq(la_feqrel(7.1824, 7.1824), DBL_MANT_DIG);
	ck_assert_int_eq(la_feqrel(INFINITY, INFINITY), DBL_MANT_DIG);
}
END_TEST


/* a few bits away from exact equality */
START_TEST (test_feqrel_few_bits_away)
{
	double w = 1;
	int i;

	for (i = 1; i < DBL_MANT_DIG - 1; ++i) {
		ck_assert_int_eq(la_feqrel(1 + w * DBL_EPSILON, 1.0),
			             DBL_MANT_DIG - i);

		ck_assert_int_eq(la_feqrel(1 - w * DBL_EPSILON, 1.0L),
			             DBL_MANT_DIG - i);

		ck_assert_int_eq(la_feqrel(1.0, 1 + (w - 1) * DBL_EPSILON),
			             DBL_MANT_DIG - i + 1);
		w *= 2;
	}

	ck_assert_int_eq(la_feqrel(1.5 + DBL_EPSILON, 1.5),
		             DBL_MANT_DIG - 1);

	ck_assert_int_eq(la_feqrel(1.5 - DBL_EPSILON, 1.5),
		             DBL_MANT_DIG - 1);

	ck_assert_int_eq(la_feqrel(1.5 - DBL_EPSILON, 1.5 + DBL_EPSILON),
		             DBL_MANT_DIG - 2);

	ck_assert_int_eq(la_feqrel(DBL_MIN / 8, DBL_MIN / 17), 3);
}
END_TEST


/* Numbers that are close */
START_TEST (test_feqrel_close_numbers)
{
	GDoubleIEEE754 x1, x2, y1, y2;

	/* 0x1.Bp+84 */
	x1.mpn.sign = 0;
	x1.mpn.mantissa_high = 0xB0000;
	x1.mpn.mantissa_low = 0;
	x1.mpn.biased_exponent = 84 + G_IEEE754_DOUBLE_BIAS;

	/* 0x1.B8p+84 */
	y1.mpn.sign = 0;
	y1.mpn.mantissa_high = 0xB8000;
	y1.mpn.mantissa_low = 0;
	y1.mpn.biased_exponent = 84 + G_IEEE754_DOUBLE_BIAS;

	ck_assert_int_eq(la_feqrel(x1.v_double, y1.v_double), 5);


	/* 0x1.8p+10 */
	x2.mpn.sign = 0;
	x2.mpn.mantissa_high = 0x80000;
	x2.mpn.mantissa_low = 0;
	x2.mpn.biased_exponent = 10 + G_IEEE754_DOUBLE_BIAS;

	/* 0x1.Cp+10 */
	y2.mpn.sign = 0;
	y2.mpn.mantissa_high = 0xC0000;
	y2.mpn.mantissa_low = 0;
	y2.mpn.biased_exponent = 10 + G_IEEE754_DOUBLE_BIAS;

	ck_assert_int_eq(la_feqrel(x2.v_double, y2.v_double), 2);


	ck_assert_int_eq(la_feqrel(1.5 * (1 - DBL_EPSILON), 1.0), 2);
	ck_assert_int_eq(la_feqrel(1.5, 1.0), 1);
	ck_assert_int_eq(la_feqrel(2 * (1 - DBL_EPSILON), 1.0), 1);
}
END_TEST

/* Factors of 2 */
START_TEST (test_feqrel_factors_of_2)
{
	ck_assert_int_eq(la_feqrel(DBL_MAX, INFINITY), 0);
	ck_assert_int_eq(la_feqrel(2 * (1 - DBL_EPSILON), 1.0), 1);
	ck_assert_int_eq(la_feqrel(1.0, 2.0), 0);
	ck_assert_int_eq(la_feqrel(4.0, 1.0), 0);
}
END_TEST

/* Extreme inequality */
START_TEST (test_feqrel_extreme_inequality)
{
	ck_assert_int_eq(la_feqrel(NAN, NAN), 0);
	ck_assert_int_eq(la_feqrel(0.0, -NAN), 0);
	ck_assert_int_eq(la_feqrel(NAN, INFINITY), 0);
	ck_assert_int_eq(la_feqrel(INFINITY, -INFINITY), 0);
	ck_assert_int_eq(la_feqrel(-DBL_MAX, INFINITY), 0);
	ck_assert_int_eq(la_feqrel(DBL_MAX, -DBL_MAX), 0);
}
END_TEST


TCase *
feqrel_tcase ()
{
    TCase *tc = tcase_create("feqrel");
    tcase_add_test(tc, test_feqrel_exact);
    tcase_add_test(tc, test_feqrel_few_bits_away);
    tcase_add_test(tc, test_feqrel_close_numbers);
    tcase_add_test(tc, test_feqrel_factors_of_2);
    tcase_add_test(tc, test_feqrel_extreme_inequality);
    return tc;
}


START_TEST (test_mknan)
{
    ck_assert(isnan(la_mknan(1)));
    ck_assert(isnan(la_mknan(LA_MAX_NAN_PAYLOAD)));
}
END_TEST

TCase *
nan_tcase ()
{
    TCase *tc = tcase_create("nan");
    tcase_add_test(tc, test_mknan);
    return tc;
}

Suite *
ieee754_suite ()
{
    Suite *s = suite_create("ieee754");
    suite_add_tcase(s, nextup_tcase());
    suite_add_tcase(s, nextdown_tcase());
    suite_add_tcase(s, ieeemean_tcase());
    suite_add_tcase(s, feqrel_tcase());
    suite_add_tcase(s, nan_tcase());
    return s;
}
