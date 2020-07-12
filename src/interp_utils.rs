pub(crate) fn bsearch(xarr: &[f64], x: f64, idx_low: usize, idx_high: usize) -> usize {
    let mut ilow = idx_low;
    let mut ihigh = idx_high;

    while ihigh > ilow + 1 {
        let i = (ihigh + ilow) / 2;
        if xarr[i] > x {
            ihigh = i;
        } else {
            ilow = i;
        }
    }
    ilow
}

pub(crate) fn linear_inter_eval(x: f64, xarr: &[f64], yarr: &[f64]) -> f64 {
    let idx = bsearch(xarr, x, 0, xarr.len() - 1);

    let x_l = xarr[idx];
    let x_h = xarr[idx + 1];
    let y_l = yarr[idx];
    let y_h = yarr[idx + 1];
    let dx = x_h - x_l;

    if dx > 0.0 {
        y_l + (x - x_l) / dx * (y_h - y_l)
    } else {
        f64::NAN
    }
}
