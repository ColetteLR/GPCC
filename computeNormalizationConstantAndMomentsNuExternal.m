function [Z, m1, m2] = computeNormalizationConstantAndMomentsNuExternal(X, Y, ...
    ma, va, mb, vb, nSplits)
    
sqrtva = sqrt(va);
	
x1 = ma - 3 * sqrtva;
x2 = ma - 1.5 * sqrtva;
x3 = ma;
x4 = ma + 1.5 * sqrtva;
x5 = ma + 3 * sqrtva;

w1 = (x5 - x1) / 90 * 7;
w2 = (x5 - x1) / 90 * 32;
w3 = (x5 - x1) / 90 * 12;
w4 = (x5 - x1) / 90 * 32;
w5 = (x5 - x1) / 90 * 7;

[z1, m11, m21] = computeNormalizationConstantAndMomentsNuInternal(X, Y, x1, mb, vb, nSplits);
[z2, m12, m22] = computeNormalizationConstantAndMomentsNuInternal(X, Y, x2, mb, vb, nSplits);
[z3, m13, m23] = computeNormalizationConstantAndMomentsNuInternal(X, Y, x3, mb, vb, nSplits);
[z4, m14, m24] = computeNormalizationConstantAndMomentsNuInternal(X, Y, x4, mb, vb, nSplits);
[z5, m15, m25] = computeNormalizationConstantAndMomentsNuInternal(X, Y, x5, mb, vb, nSplits);

z1 = z1 * (w1 * normpdf(x1, ma, sqrtva));
z2 = z2 * (w2 * normpdf(x2, ma, sqrtva));
z3 = z3 * (w3 * normpdf(x3, ma, sqrtva));
z4 = z4 * (w4 * normpdf(x4, ma, sqrtva));
z5 = z5 * (w5 * normpdf(x5, ma, sqrtva));

m11 = m11 * (w1 * normpdf(x1, ma, sqrtva));
m12 = m12 * (w2 * normpdf(x2, ma, sqrtva));
m13 = m13 * (w3 * normpdf(x3, ma, sqrtva));
m14 = m14 * (w4 * normpdf(x4, ma, sqrtva));
m15 = m15 * (w5 * normpdf(x5, ma, sqrtva));

m21 = m21 * (w1 * normpdf(x1, ma, sqrtva));
m22 = m22 * (w2 * normpdf(x2, ma, sqrtva));
m23 = m23 * (w3 * normpdf(x3, ma, sqrtva));
m24 = m24 * (w4 * normpdf(x4, ma, sqrtva));
m25 = m25 * (w5 * normpdf(x5, ma, sqrtva));

Z = z1 + z2 + z3 + z4 + z5;
m1 = m11 + m12 + m13 + m14 + m15;
m2 = m21 + m22 + m23 + m24 + m25;

end