
#include <math.h>
#include <float.h>
#include <glib.h>
#include <check.h>
#include <libla/ieee754.h>

/* eqrel tests ported from Tango v 0.99.9,
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

Suite *
ieee754_suite ()
{
    Suite *s = suite_create("ieee754");
    suite_add_tcase(s, feqrel_tcase());
    return s;
}
