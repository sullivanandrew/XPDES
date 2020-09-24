#include <math.h>

double KERRfISOout (
  double r,
  double theta,
  double r_H,
  double M)
{
  return(pow(0.1e1 - 0.1e1 * r_H * r_H * pow(r, -0.2e1), 0.2e1) * (0.2e1 * M * M * pow(r, -0.2e1) + pow(0.1e1 - 0.1e1 * r_H * r_H * pow(r, -0.2e1), 0.2e1) + 0.2e1 * M / r * (0.1e1 + r_H * r_H * pow(r, -0.2e1)) - 0.1e1 * (M * M - 0.4e1 * r_H * r_H) * pow(r, -0.2e1) * pow(sin(theta), 0.2e1)) / (pow(0.2e1 * M * M * pow(r, -0.2e1) + pow(0.1e1 - 0.1e1 * r_H * r_H * pow(r, -0.2e1), 0.2e1) + 0.2e1 * M / r * (0.1e1 + r_H * r_H * pow(r, -0.2e1)), 0.2e1) - 0.1e1 * pow(0.1e1 - 0.1e1 * r_H * r_H * pow(r, -0.2e1), 0.2e1) * (M * M - 0.4e1 * r_H * r_H) * pow(r, -0.2e1) * pow(sin(theta), 0.2e1)));
}
#include <math.h>

double KERRmISOout (
  double r,
  double theta,
  double r_H,
  double M)
{
  return(pow(0.1e1 - 0.1e1 * r_H * r_H * pow(r, -0.2e1), 0.2e1) * pow(0.2e1 * M * M * pow(r, -0.2e1) + pow(0.1e1 - 0.1e1 * r_H * r_H * pow(r, -0.2e1), 0.2e1) + 0.2e1 * M / r * (0.1e1 + r_H * r_H * pow(r, -0.2e1)) - 0.1e1 * (M * M - 0.4e1 * r_H * r_H) * pow(r, -0.2e1) * pow(sin(theta), 0.2e1), 0.2e1) / (pow(0.2e1 * M * M * pow(r, -0.2e1) + pow(0.1e1 - 0.1e1 * r_H * r_H * pow(r, -0.2e1), 0.2e1) + 0.2e1 * M / r * (0.1e1 + r_H * r_H * pow(r, -0.2e1)), 0.2e1) - 0.1e1 * pow(0.1e1 - 0.1e1 * r_H * r_H * pow(r, -0.2e1), 0.2e1) * (M * M - 0.4e1 * r_H * r_H) * pow(r, -0.2e1) * pow(sin(theta), 0.2e1)));
}
#include <math.h>

double KERRlISOout (
  double r,
  double theta,
  double r_H,
  double M)
{
  return(pow(0.1e1 - 0.1e1 * r_H * r_H * pow(r, -0.2e1), 0.2e1));
}
#include <math.h>

double KERRomegaISOout (
  double r,
  double theta,
  double r_H,
  double M)
{
  return(0.2e1 * M * sqrt(M * M - 0.4e1 * r_H * r_H) * pow(r, -0.2e1) * (0.1e1 + M / r + r_H * r_H * pow(r, -0.2e1)) / (pow(0.2e1 * M * M * pow(r, -0.2e1) + pow(0.1e1 - 0.1e1 * r_H * r_H * pow(r, -0.2e1), 0.2e1) + 0.2e1 * M / r * (0.1e1 + r_H * r_H * pow(r, -0.2e1)), 0.2e1) - 0.1e1 * pow(0.1e1 - 0.1e1 * r_H * r_H * pow(r, -0.2e1), 0.2e1) * (M * M - 0.4e1 * r_H * r_H) * pow(r, -0.2e1) * pow(sin(theta), 0.2e1)));
}
#include <math.h>

double KERRu0Xout (
  double x,
  double y,
  double r_H,
  double M)
{
  return(-0.1e1 * pow(x - 0.2e1, 0.2e1) * r_H * r_H * (pow(x - 0.1e1, 0.2e1) * (M - 0.2e1 * r_H) * (M + 0.2e1 * r_H) * pow(sin(y), 0.2e1) - 0.1e1 * pow(x, 0.4e1) * r_H * r_H + 0.2e1 * r_H * (M + 0.2e1 * r_H) * pow(x, 0.3e1) - 0.2e1 * (M + 0.2e1 * r_H) * (M + r_H) * x * x + 0.4e1 * M * (M + 0.2e1 * r_H) * x - 0.2e1 * M * (M + 0.2e1 * r_H)) * x * x / (-0.1e1 * x * x * r_H * r_H * pow(x - 0.1e1, 0.2e1) * pow(x - 0.2e1, 0.2e1) * (M - 0.2e1 * r_H) * (M + 0.2e1 * r_H) * pow(sin(y), 0.2e1) + 0.4e1 * pow(0.5000000000e0 * pow(x, 0.4e1) * r_H * r_H - 0.1e1 * r_H * (M + 0.2e1 * r_H) * pow(x, 0.3e1) + (M + 0.2e1 * r_H) * (M + r_H) * x * x - 0.2e1 * M * (M + 0.2e1 * r_H) * x + M * (M + 0.2e1 * r_H), 0.2e1)));
}
#include <math.h>

double KERRu1Xout (
  double x,
  double y,
  double r_H,
  double M)
{
  return(pow(x - 0.2e1, 0.2e1) * pow(pow(x - 0.1e1, 0.2e1) * (M - 0.2e1 * r_H) * (M + 0.2e1 * r_H) * pow(sin(y), 0.2e1) - 0.1e1 * pow(x, 0.4e1) * r_H * r_H + 0.2e1 * r_H * (M + 0.2e1 * r_H) * pow(x, 0.3e1) - 0.2e1 * (M + 0.2e1 * r_H) * (M + r_H) * x * x + 0.4e1 * M * (M + 0.2e1 * r_H) * x - 0.2e1 * M * (M + 0.2e1 * r_H), 0.2e1) * x * x / (-0.1e1 * x * x * r_H * r_H * pow(x - 0.1e1, 0.2e1) * pow(x - 0.2e1, 0.2e1) * (M - 0.2e1 * r_H) * (M + 0.2e1 * r_H) * pow(sin(y), 0.2e1) + 0.4e1 * pow(0.5000000000e0 * pow(x, 0.4e1) * r_H * r_H - 0.1e1 * r_H * (M + 0.2e1 * r_H) * pow(x, 0.3e1) + (M + 0.2e1 * r_H) * (M + r_H) * x * x - 0.2e1 * M * (M + 0.2e1 * r_H) * x + M * (M + 0.2e1 * r_H), 0.2e1)));
}
#include <math.h>

double KERRu2Xout (
  double x,
  double y,
  double r_H,
  double M)
{
  return(x * x * pow(x - 0.2e1, 0.2e1));
}
#include <math.h>

double KERRu3Xout (
  double x,
  double y,
  double r_H,
  double M)
{
  return(-0.2e1 * pow(x - 0.1e1, 0.2e1) * (-0.1e1 * x * x * r_H + (M + 0.2e1 * r_H) * x - 0.1e1 * M - 0.2e1 * r_H) * sqrt(M * M - 0.4e1 * r_H * r_H) * r_H * M / (-0.1e1 * x * x * r_H * r_H * pow(x - 0.1e1, 0.2e1) * pow(x - 0.2e1, 0.2e1) * (M - 0.2e1 * r_H) * (M + 0.2e1 * r_H) * pow(sin(y), 0.2e1) + 0.4e1 * pow(0.5000000000e0 * pow(x, 0.4e1) * r_H * r_H - 0.1e1 * r_H * (M + 0.2e1 * r_H) * pow(x, 0.3e1) + (M + 0.2e1 * r_H) * (M + r_H) * x * x - 0.2e1 * M * (M + 0.2e1 * r_H) * x + M * (M + 0.2e1 * r_H), 0.2e1)));
}
double KERRu4Xout (
  double x,
  double y,
  double r_H,
  double M)
{
  return(0.0e0);
}
#include <math.h>

double KERRu5Xout (
  double x,
  double y,
  double r_H,
  double M)
{
  return(-0.5000000000e0 * cos(y) * pow(x - 0.1e1, 0.3e1) * (-0.1e1 * x * x * r_H + (M + 0.2e1 * r_H) * x - 0.1e1 * M - 0.2e1 * r_H) * pow(x - 0.2e1, 0.2e1) * sin(y) * (M + 0.2e1 * r_H) * (M - 0.2e1 * r_H) * r_H * r_H * M * (0.5000000000e0 * pow(x, 0.4e1) * r_H * r_H - 0.1e1 * r_H * (M + 0.2e1 * r_H) * pow(x, 0.3e1) + (M + 0.2e1 * r_H) * (M + r_H) * x * x - 0.2e1 * M * (M + 0.2e1 * r_H) * x + M * (M + 0.2e1 * r_H)) * x * x * pow(-0.2500000000e0 * x * x * r_H * r_H * pow(x - 0.1e1, 0.2e1) * pow(x - 0.2e1, 0.2e1) * (M - 0.2e1 * r_H) * (M + 0.2e1 * r_H) * pow(sin(y), 0.2e1) + pow(0.5000000000e0 * pow(x, 0.4e1) * r_H * r_H - 0.1e1 * r_H * (M + 0.2e1 * r_H) * pow(x, 0.3e1) + (M + 0.2e1 * r_H) * (M + r_H) * x * x - 0.2e1 * M * (M + 0.2e1 * r_H) * x + M * (M + 0.2e1 * r_H), 0.2e1), -0.2e1));
}
#include <math.h>

double KERRu6Xout (
  double x,
  double y,
  double r_H,
  double M)
{
  return(cos(y) * pow(x - 0.1e1, 0.2e1) * pow(x - 0.2e1, 0.2e1) * sin(y) * (M + 0.2e1 * r_H) * (pow(x - 0.1e1, 0.2e1) * (M - 0.2e1 * r_H) * (M + 0.2e1 * r_H) * pow(sin(y), 0.2e1) - 0.1e1 * pow(x, 0.4e1) * r_H * r_H + 0.2e1 * r_H * (M + 0.2e1 * r_H) * pow(x, 0.3e1) - 0.2e1 * (M + 0.2e1 * r_H) * (M + r_H) * x * x + 0.4e1 * M * (M + 0.2e1 * r_H) * x - 0.2e1 * M * (M + 0.2e1 * r_H)) * (M - 0.2e1 * r_H) * (-0.1250000000e0 * x * x * r_H * r_H * pow(x - 0.1e1, 0.2e1) * pow(x - 0.2e1, 0.2e1) * (M - 0.2e1 * r_H) * (M + 0.2e1 * r_H) * pow(sin(y), 0.2e1) + (-0.5000000000e0 * x * x * r_H + (M + 0.2e1 * r_H) * x - 0.1e1 * M - 0.2e1 * r_H) * (-0.5000000000e0 * x * x * r_H + M * x - 0.1e1 * M) * (0.5000000000e0 * pow(x, 0.4e1) * r_H * r_H - 0.1e1 * r_H * (M + 0.2e1 * r_H) * pow(x, 0.3e1) + (M + 0.2e1 * r_H) * (M + r_H) * x * x - 0.2e1 * M * (M + 0.2e1 * r_H) * x + M * (M + 0.2e1 * r_H))) * x * x * pow(-0.2500000000e0 * x * x * r_H * r_H * pow(x - 0.1e1, 0.2e1) * pow(x - 0.2e1, 0.2e1) * (M - 0.2e1 * r_H) * (M + 0.2e1 * r_H) * pow(sin(y), 0.2e1) + pow(0.5000000000e0 * pow(x, 0.4e1) * r_H * r_H - 0.1e1 * r_H * (M + 0.2e1 * r_H) * pow(x, 0.3e1) + (M + 0.2e1 * r_H) * (M + r_H) * x * x - 0.2e1 * M * (M + 0.2e1 * r_H) * x + M * (M + 0.2e1 * r_H), 0.2e1), -0.2e1));
}
double KERRu7Xout (
  double x,
  double y,
  double r_H,
  double M)
{
  return(0.0e0);
}
#include <math.h>

double KERRu8Xout (
  double x,
  double y,
  double r_H,
  double M)
{
  return(-0.2500000000e0 * cos(y) * pow(x - 0.1e1, 0.4e1) * (-0.1e1 * x * x * r_H + (M + 0.2e1 * r_H) * x - 0.1e1 * M - 0.2e1 * r_H) * pow(x - 0.2e1, 0.2e1) * sin(y) * pow(r_H, 0.3e1) * M * pow(M * M - 0.4e1 * r_H * r_H, 0.3e1 / 0.2e1) * x * x * pow(-0.2500000000e0 * x * x * r_H * r_H * pow(x - 0.1e1, 0.2e1) * pow(x - 0.2e1, 0.2e1) * (M - 0.2e1 * r_H) * (M + 0.2e1 * r_H) * pow(sin(y), 0.2e1) + pow(0.5000000000e0 * pow(x, 0.4e1) * r_H * r_H - 0.1e1 * r_H * (M + 0.2e1 * r_H) * pow(x, 0.3e1) + (M + 0.2e1 * r_H) * (M + r_H) * x * x - 0.2e1 * M * (M + 0.2e1 * r_H) * x + M * (M + 0.2e1 * r_H), 0.2e1), -0.2e1));
}
double KERRu9Xout (
  double x,
  double y,
  double r_H,
  double M)
{
  return(0.0e0);
}
#include <math.h>

double KERRdu0dxXout (
  double x,
  double y,
  double r_H,
  double M)
{
  return(0.5000000000e0 * r_H * r_H * ((M + 0.2e1 * r_H) * pow(x - 0.1e1, 0.2e1) * (M - 0.2e1 * r_H) * (0.2500000000e0 * pow(r_H, 0.3e1) * pow(x, 0.8e1) - 0.2e1 * pow(r_H, 0.3e1) * pow(x, 0.7e1) + (-0.1e1 * M * M * r_H + 0.6e1 * pow(r_H, 0.3e1)) * pow(x, 0.6e1) + (M + 0.2e1 * r_H) * (M * M + 0.4e1 * M * r_H - 0.4e1 * r_H * r_H) * pow(x, 0.5e1) + (-0.5e1 * pow(M, 0.3e1) - 0.21e2 * M * M * r_H - 0.20e2 * M * r_H * r_H + 0.4e1 * pow(r_H, 0.3e1)) * pow(x, 0.4e1) + 0.11e2 * M * pow(M + 0.2e1 * r_H, 0.2e1) * pow(x, 0.3e1) - 0.13e2 * M * pow(M + 0.2e1 * r_H, 0.2e1) * x * x + 0.8e1 * M * pow(M + 0.2e1 * r_H, 0.2e1) * x - 0.2e1 * M * pow(M + 0.2e1 * r_H, 0.2e1)) * pow(cos(y), 0.2e1) + (-0.2500000000e0 * pow(r_H, 0.3e1) * pow(x, 0.8e1) + r_H * r_H * (M + 0.2e1 * r_H) * pow(x, 0.7e1) - 0.2e1 * (M + 0.2e1 * r_H) * (M + 0.1500000000e1 * r_H) * r_H * pow(x, 0.6e1) + (M + 0.2e1 * r_H) * (M * M + 0.10e2 * M * r_H + 0.4e1 * r_H * r_H) * pow(x, 0.5e1) + (-0.5e1 * pow(M, 0.3e1) - 0.34e2 * M * M * r_H - 0.50e2 * M * r_H * r_H - 0.4e1 * pow(r_H, 0.3e1)) * pow(x, 0.4e1) + 0.11e2 * (M + 0.2e1 * r_H) * (M + 0.3090909091e1 * r_H) * M * pow(x, 0.3e1) - 0.13e2 * (M + 0.2e1 * r_H) * M * (M + 0.2307692308e1 * r_H) * x * x + 0.8e1 * M * pow(M + 0.2e1 * r_H, 0.2e1) * x - 0.2e1 * M * pow(M + 0.2e1 * r_H, 0.2e1)) * pow(-0.1e1 * x * x * r_H + (M + 0.2e1 * r_H) * x - 0.1e1 * M - 0.2e1 * r_H, 0.2e1)) * M * (x - 0.2e1) * x * pow(0.2500000000e0 * x * x * r_H * r_H * pow(x - 0.1e1, 0.2e1) * pow(x - 0.2e1, 0.2e1) * (M - 0.2e1 * r_H) * (M + 0.2e1 * r_H) * pow(cos(y), 0.2e1) + (-0.2500000000e0 * pow(x, 0.6e1) * pow(r_H, 0.3e1) + 0.7500000000e0 * r_H * r_H * (M + 0.2e1 * r_H) * pow(x, 0.5e1) - 0.1e1 * (M + 0.2e1 * r_H) * (M + 0.1750000000e1 * r_H) * r_H * pow(x, 0.4e1) + (M + 0.2e1 * r_H) * (M * M + 0.2e1 * M * r_H + 0.2e1 * r_H * r_H) * pow(x, 0.3e1) + (-0.3e1 * pow(M, 0.3e1) - 0.7e1 * M * M * r_H - 0.3e1 * M * r_H * r_H - 0.2e1 * pow(r_H, 0.3e1)) * x * x + 0.3e1 * M * M * (M + 0.2e1 * r_H) * x - 0.1e1 * M * M * (M + 0.2e1 * r_H)) * (-0.1e1 * x * x * r_H + (M + 0.2e1 * r_H) * x - 0.1e1 * M - 0.2e1 * r_H), -0.2e1));
}
#include <math.h>

double KERRdu1dxXout (
  double x,
  double y,
  double r_H,
  double M)
{
  return((x - 0.1e1) * (pow(x - 0.1e1, 0.2e1) * (M - 0.2e1 * r_H) * (M + 0.2e1 * r_H) * pow(cos(y), 0.2e1) + pow(-0.1e1 * x * x * r_H + (M + 0.2e1 * r_H) * x - 0.1e1 * M - 0.2e1 * r_H, 0.2e1)) * (0.1250000000e0 * r_H * r_H * pow(x, 0.3e1) * pow(x - 0.1e1, 0.2e1) * pow(x - 0.2e1, 0.3e1) * pow(M - 0.2e1 * r_H, 0.2e1) * pow(M + 0.2e1 * r_H, 0.2e1) * pow(cos(y), 0.4e1) + (M + 0.2e1 * r_H) * (M - 0.2e1 * r_H) * (0.3750000000e0 * pow(x, 0.10e2) * pow(r_H, 0.4e1) - 0.7500000000e0 * pow(r_H, 0.3e1) * (M + 0.5e1 * r_H) * pow(x, 0.9e1) + r_H * r_H * (0.1625000000e2 * r_H * r_H + M * M + 0.6750000000e1 * M * r_H) * pow(x, 0.8e1) + (-0.40e2 * pow(r_H, 0.4e1) - 0.2550000000e2 * pow(r_H, 0.3e1) * M - 0.1500000000e1 * pow(M, 0.3e1) * r_H - 0.8e1 * M * M * r_H * r_H) * pow(x, 0.7e1) + (0.61e2 * pow(r_H, 0.4e1) + 0.29e2 * M * M * r_H * r_H + 0.5250000000e2 * pow(r_H, 0.3e1) * M + 0.1050000000e2 * pow(M, 0.3e1) * r_H + pow(M, 0.4e1)) * pow(x, 0.6e1) + (-0.58e2 * pow(r_H, 0.4e1) - 0.3450000000e2 * pow(M, 0.3e1) * r_H - 0.62e2 * M * M * r_H * r_H - 0.63e2 * pow(r_H, 0.3e1) * M - 0.6e1 * pow(M, 0.4e1)) * pow(x, 0.5e1) + (0.6750000000e2 * pow(M, 0.3e1) * r_H + 0.88e2 * M * M * r_H * r_H + 0.42e2 * pow(r_H, 0.3e1) * M + 0.32e2 * pow(r_H, 0.4e1) + 0.15e2 * pow(M, 0.4e1)) * pow(x, 0.4e1) + (-0.20e2 * pow(M, 0.4e1) - 0.82e2 * pow(M, 0.3e1) * r_H - 0.88e2 * M * M * r_H * r_H - 0.12e2 * pow(r_H, 0.3e1) * M - 0.8e1 * pow(r_H, 0.4e1)) * pow(x, 0.3e1) + 0.15e2 * M * M * pow(M + 0.2e1 * r_H, 0.2e1) * x * x - 0.6e1 * M * M * pow(M + 0.2e1 * r_H, 0.2e1) * x + M * M * pow(M + 0.2e1 * r_H, 0.2e1)) * pow(cos(y), 0.2e1) + pow(-0.1e1 * x * x * r_H + (M + 0.2e1 * r_H) * x - 0.1e1 * M - 0.2e1 * r_H, 0.2e1) * (0.2500000000e0 * pow(r_H, 0.4e1) * pow(x, 0.8e1) - 0.1e1 * pow(r_H, 0.3e1) * (M + 0.2e1 * r_H) * pow(x, 0.7e1) + 0.1875000000e1 * (M + 0.2e1 * r_H) * (M + 0.1733333333e1 * r_H) * r_H * r_H * pow(x, 0.6e1) + (-0.11e2 * pow(r_H, 0.4e1) - 0.2500000000e1 * pow(M, 0.3e1) * r_H - 0.1125000000e2 * M * M * r_H * r_H - 0.18e2 * pow(r_H, 0.3e1) * M) * pow(x, 0.5e1) + (pow(M, 0.3e1) + 0.1050000000e2 * M * M * r_H + 0.7500000000e1 * M * r_H * r_H + 0.5e1 * pow(r_H, 0.3e1)) * (M + 0.2e1 * r_H) * pow(x, 0.4e1) + (-0.4e1 * pow(M, 0.4e1) - 0.26e2 * pow(M, 0.3e1) * r_H - 0.39e2 * M * M * r_H * r_H - 0.8e1 * pow(r_H, 0.3e1) * M - 0.4e1 * pow(r_H, 0.4e1)) * pow(x, 0.3e1) + 0.6e1 * (M + 0.2e1 * r_H) * (M + 0.2666666667e1 * r_H) * M * M * x * x - 0.4e1 * M * M * pow(M + 0.2e1 * r_H, 0.2e1) * x + M * M * pow(M + 0.2e1 * r_H, 0.2e1))) * (x - 0.2e1) * x * pow(0.2500000000e0 * x * x * r_H * r_H * pow(x - 0.1e1, 0.2e1) * pow(x - 0.2e1, 0.2e1) * (M - 0.2e1 * r_H) * (M + 0.2e1 * r_H) * pow(cos(y), 0.2e1) + (-0.2500000000e0 * pow(x, 0.6e1) * pow(r_H, 0.3e1) + 0.7500000000e0 * r_H * r_H * (M + 0.2e1 * r_H) * pow(x, 0.5e1) - 0.1e1 * (M + 0.2e1 * r_H) * (M + 0.1750000000e1 * r_H) * r_H * pow(x, 0.4e1) + (M + 0.2e1 * r_H) * (M * M + 0.2e1 * M * r_H + 0.2e1 * r_H * r_H) * pow(x, 0.3e1) + (-0.3e1 * pow(M, 0.3e1) - 0.7e1 * M * M * r_H - 0.3e1 * M * r_H * r_H - 0.2e1 * pow(r_H, 0.3e1)) * x * x + 0.3e1 * M * M * (M + 0.2e1 * r_H) * x - 0.1e1 * M * M * (M + 0.2e1 * r_H)) * (-0.1e1 * x * x * r_H + (M + 0.2e1 * r_H) * x - 0.1e1 * M - 0.2e1 * r_H), -0.2e1));
}
#include <math.h>

double KERRdu2dxXout (
  double x,
  double y,
  double r_H,
  double M)
{
  return(0.4e1 * pow(x, 0.3e1) - 0.12e2 * x * x + 0.8e1 * x);
}
#include <math.h>

double KERRdu3dxXout (
  double x,
  double y,
  double r_H,
  double M)
{
  return(0.5000000000e0 * (x - 0.1e1) * (0.7500000000e0 * (M + 0.2e1 * r_H) * pow(x - 0.1e1, 0.3e1) * (M - 0.2e1 * r_H) * (-0.6666666667e0 * r_H * pow(x, 0.3e1) + (M + 0.2e1 * r_H) * x * x + (-0.2e1 * M - 0.4e1 * r_H) * x + 0.1333333333e1 * M + 0.2666666667e1 * r_H) * r_H * r_H * (x - 0.2e1) * x * pow(cos(y), 0.2e1) + pow(-0.1e1 * x * x * r_H + (M + 0.2e1 * r_H) * x - 0.1e1 * M - 0.2e1 * r_H, 0.2e1) * (-0.1e1 * pow(x, 0.6e1) * pow(r_H, 0.3e1) + (0.2250000000e1 * M * r_H * r_H + 0.6e1 * pow(r_H, 0.3e1)) * pow(x, 0.5e1) - 0.2e1 * (M + 0.2e1 * r_H) * (M + 0.3625000000e1 * r_H) * r_H * pow(x, 0.4e1) + (M + 0.2e1 * r_H) * pow(M + 0.3e1 * r_H, 0.2e1) * pow(x, 0.3e1) - 0.3e1 * (M + 0.2e1 * r_H) * (M * M + 0.2e1 * M * r_H + 0.2e1 * r_H * r_H) * x * x + (0.3e1 * pow(M, 0.3e1) + 0.8e1 * M * M * r_H + 0.6e1 * M * r_H * r_H + 0.4e1 * pow(r_H, 0.3e1)) * x - 0.1e1 * M * M * (M + 0.2e1 * r_H))) * sqrt(M * M - 0.4e1 * r_H * r_H) * r_H * M * pow(0.2500000000e0 * x * x * r_H * r_H * pow(x - 0.1e1, 0.2e1) * pow(x - 0.2e1, 0.2e1) * (M - 0.2e1 * r_H) * (M + 0.2e1 * r_H) * pow(cos(y), 0.2e1) + (-0.2500000000e0 * pow(x, 0.6e1) * pow(r_H, 0.3e1) + 0.7500000000e0 * r_H * r_H * (M + 0.2e1 * r_H) * pow(x, 0.5e1) - 0.1e1 * (M + 0.2e1 * r_H) * (M + 0.1750000000e1 * r_H) * r_H * pow(x, 0.4e1) + (M + 0.2e1 * r_H) * (M * M + 0.2e1 * M * r_H + 0.2e1 * r_H * r_H) * pow(x, 0.3e1) + (-0.3e1 * pow(M, 0.3e1) - 0.7e1 * M * M * r_H - 0.3e1 * M * r_H * r_H - 0.2e1 * pow(r_H, 0.3e1)) * x * x + 0.3e1 * M * M * (M + 0.2e1 * r_H) * x - 0.1e1 * M * M * (M + 0.2e1 * r_H)) * (-0.1e1 * x * x * r_H + (M + 0.2e1 * r_H) * x - 0.1e1 * M - 0.2e1 * r_H), -0.2e1));
}
double KERRdu4dxXout (
  double x,
  double y,
  double r_H,
  double M)
{
  return(0.0e0);
}
#include <math.h>

double KERRd2u0dx2Xout (
  double x,
  double y,
  double r_H,
  double M)
{
  return(0.5000000000e0 * (-0.7500000000e0 * pow(M + 0.2e1 * r_H, 0.2e1) * pow(x - 0.1e1, 0.3e1) * pow(M - 0.2e1 * r_H, 0.2e1) * r_H * r_H * ((-0.1666666667e0 * pow(r_H, 0.3e1) - 0.6666666667e0 * M * M * r_H) * pow(x, 0.8e1) + (0.1333333333e1 * pow(r_H, 0.3e1) + pow(M, 0.3e1) + 0.5333333333e1 * M * M * r_H + 0.4e1 * M * r_H * r_H) * pow(x, 0.7e1) + (-0.7e1 * pow(M, 0.3e1) - 0.28e2 * M * M * r_H - 0.28e2 * M * r_H * r_H - 0.4e1 * pow(r_H, 0.3e1)) * pow(x, 0.6e1) + 0.2366666667e2 * (M + 0.2e1 * r_H) * (M * M + 0.1943661972e1 * M * r_H + 0.1126760563e0 * r_H * r_H) * pow(x, 0.5e1) - 0.4833333333e2 * (M + 0.2e1 * r_H) * (M * M + 0.1986206897e1 * M * r_H + 0.2758620690e-1 * r_H * r_H) * pow(x, 0.4e1) + 0.62e2 * M * pow(M + 0.2e1 * r_H, 0.2e1) * pow(x, 0.3e1) - 0.4866666667e2 * M * pow(M + 0.2e1 * r_H, 0.2e1) * x * x + 0.2133333333e2 * M * pow(M + 0.2e1 * r_H, 0.2e1) * x - 0.4e1 * M * pow(M + 0.2e1 * r_H, 0.2e1)) * pow(x - 0.2e1, 0.2e1) * x * x * pow(cos(y), 0.4e1) + (-0.3750000000e0 * pow(x, 0.18e2) * pow(r_H, 0.7e1) + 0.8750000000e0 * (0.7714285714e1 * r_H + M) * pow(r_H, 0.6e1) * pow(x, 0.17e2) + (-0.5475000000e2 * pow(r_H, 0.7e1) + 0.1500000000e1 * M * M * pow(r_H, 0.5e1) - 0.1487500000e2 * M * pow(r_H, 0.6e1)) * pow(x, 0.16e2) + (0.264e3 * pow(r_H, 0.7e1) + 0.1032500000e3 * M * pow(r_H, 0.6e1) - 0.8250000000e1 * pow(M, 0.3e1) * pow(r_H, 0.4e1) - 0.24e2 * M * M * pow(r_H, 0.5e1)) * pow(x, 0.15e2) + (-0.840e3 * pow(r_H, 0.7e1) + 0.2287500000e3 * M * M * pow(r_H, 0.5e1) - 0.3587500000e3 * M * pow(r_H, 0.6e1) + 0.1150000000e2 * pow(M, 0.4e1) * pow(r_H, 0.3e1) + 0.1237500000e3 * pow(M, 0.3e1) * pow(r_H, 0.4e1)) * pow(x, 0.14e2) - 0.6e1 * (M + 0.2e1 * r_H) * (pow(M, 0.4e1) + 0.2483333333e2 * pow(M, 0.3e1) * r_H + 0.1067083333e3 * M * M * r_H * r_H + 0.4033333333e2 * pow(r_H, 0.3e1) * M - 0.154e3 * pow(r_H, 0.4e1)) * r_H * r_H * pow(x, 0.13e2) + 0.78e2 * (M + 0.2e1 * r_H) * r_H * r_H * (pow(M, 0.4e1) + 0.1201282051e2 * pow(M, 0.3e1) * r_H + 0.3609935897e2 * M * M * r_H * r_H + 0.1944551282e2 * pow(r_H, 0.3e1) * M - 0.1830769231e2 * pow(r_H, 0.4e1)) * pow(x, 0.12e2) + (0.3072e4 * pow(r_H, 0.7e1) - 0.476e3 * pow(M, 0.5e1) * r_H * r_H - 0.4744e4 * pow(M, 0.4e1) * pow(r_H, 0.3e1) - 0.1680250000e5 * pow(M, 0.3e1) * pow(r_H, 0.4e1) - 0.24081e5 * M * M * pow(r_H, 0.5e1) - 0.9816e4 * M * pow(r_H, 0.6e1) + pow(M, 0.7e1)) * pow(x, 0.11e2) + (-0.2208e4 * pow(r_H, 0.7e1) - 0.6e1 * pow(M, 0.6e1) * r_H + 0.1804e4 * pow(M, 0.5e1) * r_H * r_H + 0.1447550000e5 * pow(M, 0.4e1) * pow(r_H, 0.3e1) + 0.4454450000e5 * pow(M, 0.3e1) * pow(r_H, 0.4e1) + 0.59577e5 * M * M * pow(r_H, 0.5e1) + 0.27324e5 * M * pow(r_H, 0.6e1) - 0.11e2 * pow(M, 0.7e1)) * pow(x, 0.10e2) + (0.55e2 * pow(M, 0.7e1) + 0.60e2 * pow(M, 0.6e1) * r_H - 0.4646e4 * pow(M, 0.5e1) * r_H * r_H - 0.32203e5 * pow(M, 0.4e1) * pow(r_H, 0.3e1) - 0.88582e5 * pow(M, 0.3e1) * pow(r_H, 0.4e1) - 0.110340e6 * M * M * pow(r_H, 0.5e1) - 0.50760e5 * M * pow(r_H, 0.6e1) + 0.960e3 * pow(r_H, 0.7e1)) * pow(x, 0.9e1) + (-0.165e3 * pow(M, 0.7e1) - 0.290e3 * pow(M, 0.6e1) * r_H + 0.8286e4 * pow(M, 0.5e1) * r_H * r_H + 0.52728e5 * pow(M, 0.4e1) * pow(r_H, 0.3e1) + 0.132717e6 * pow(M, 0.3e1) * pow(r_H, 0.4e1) + 0.154708e6 * M * M * pow(r_H, 0.5e1) + 0.68980e5 * M * pow(r_H, 0.6e1) - 0.192e3 * pow(r_H, 0.7e1)) * pow(x, 0.8e1) + 0.333e3 * pow(M + 0.2e1 * r_H, 0.2e1) * (pow(M, 0.4e1) - 0.1357357357e1 * pow(M, 0.3e1) * r_H - 0.2826426426e2 * M * M * r_H * r_H - 0.7041441441e2 * pow(r_H, 0.3e1) * M - 0.5310510511e2 * pow(r_H, 0.4e1)) * M * pow(x, 0.7e1) - 0.483e3 * pow(M + 0.2e1 * r_H, 0.2e1) * (pow(M, 0.4e1) - 0.2194616977e0 * pow(M, 0.3e1) * r_H - 0.1713457557e2 * M * M * r_H * r_H - 0.3967701863e2 * pow(r_H, 0.3e1) * M - 0.2838095238e2 * pow(r_H, 0.4e1)) * M * pow(x, 0.6e1) + 0.525e3 * pow(M + 0.2e1 * r_H, 0.3e1) * (pow(M, 0.3e1) - 0.9180952381e0 * M * M * r_H - 0.7116190476e1 * M * r_H * r_H - 0.7481904762e1 * pow(r_H, 0.3e1)) * M * pow(x, 0.5e1) - 0.435e3 * pow(M + 0.2e1 * r_H, 0.3e1) * (pow(M, 0.3e1) + 0.3218390805e0 * M * M * r_H - 0.3577011494e1 * M * r_H * r_H - 0.3632183908e1 * pow(r_H, 0.3e1)) * M * pow(x, 0.4e1) + 0.270e3 * (pow(M, 0.3e1) + 0.1259259259e1 * M * M * r_H - 0.1481481481e1 * M * r_H * r_H - 0.1481481481e1 * pow(r_H, 0.3e1)) * pow(M + 0.2e1 * r_H, 0.3e1) * M * pow(x, 0.3e1) - 0.118e3 * pow(M + 0.2e1 * r_H, 0.3e1) * M * (pow(M, 0.3e1) + 0.1796610169e1 * M * M * r_H - 0.4067796610e0 * M * r_H * r_H - 0.4067796610e0 * pow(r_H, 0.3e1)) * x * x + 0.32e2 * pow(M, 0.3e1) * pow(M + 0.2e1 * r_H, 0.4e1) * x - 0.4e1 * pow(M, 0.3e1) * pow(M + 0.2e1 * r_H, 0.4e1)) * (M + 0.2e1 * r_H) * (x - 0.1e1) * (M - 0.2e1 * r_H) * pow(cos(y), 0.2e1) + pow(-0.1e1 * x * x * r_H + (M + 0.2e1 * r_H) * x - 0.1e1 * M - 0.2e1 * r_H, 0.3e1) * (-0.1250000000e0 * pow(x, 0.15e2) * pow(r_H, 0.6e1) + 0.7500000000e0 * (M + 0.2500000000e1 * r_H) * pow(r_H, 0.5e1) * pow(x, 0.14e2) - 0.2250000000e1 * (M + 0.2e1 * r_H) * (M + 0.2666666667e1 * r_H) * pow(r_H, 0.4e1) * pow(x, 0.13e2) + 0.2250000000e1 * (M * M + 0.11e2 * M * r_H + 0.9388888889e1 * r_H * r_H) * (M + 0.2e1 * r_H) * pow(r_H, 0.3e1) * pow(x, 0.12e2) + 0.7500000000e0 * r_H * r_H * (M + 0.2e1 * r_H) * (pow(M, 0.3e1) - 0.38e2 * M * M * r_H - 0.172e3 * M * r_H * r_H - 0.58e2 * pow(r_H, 0.3e1)) * pow(x, 0.11e2) - 0.3e1 * (M + 0.2e1 * r_H) * (pow(M, 0.4e1) + 0.7500000000e0 * pow(M, 0.3e1) * r_H - 0.5475000000e2 * M * M * r_H * r_H - 0.1435000000e3 * pow(r_H, 0.3e1) * M - 0.1650000000e2 * pow(r_H, 0.4e1)) * r_H * pow(x, 0.10e2) + (M + 0.2e1 * r_H) * (pow(M, 0.5e1) + 0.28e2 * pow(M, 0.4e1) * r_H - 0.1025000000e2 * pow(M, 0.3e1) * r_H * r_H - 0.587e3 * M * M * pow(r_H, 0.3e1) - 0.1017e4 * M * pow(r_H, 0.4e1) - 0.20e2 * pow(r_H, 0.5e1)) * pow(x, 0.9e1) - 0.9e1 * (M + 0.2e1 * r_H) * (pow(M, 0.5e1) + 0.1333333333e2 * pow(M, 0.4e1) * r_H - 0.8416666667e1 * pow(M, 0.3e1) * r_H * r_H - 0.1603333333e3 * M * M * pow(r_H, 0.3e1) - 0.1983333333e3 * M * pow(r_H, 0.4e1) + 0.2e1 * pow(r_H, 0.5e1)) * pow(x, 0.8e1) + 0.36e2 * (pow(M, 0.4e1) + 0.6666666667e1 * pow(M, 0.3e1) * r_H - 0.1870833333e2 * M * M * r_H * r_H - 0.3316666667e2 * pow(r_H, 0.3e1) * M + 0.3333333333e0 * pow(r_H, 0.4e1)) * pow(M + 0.2e1 * r_H, 0.2e1) * pow(x, 0.7e1) - 0.84e2 * pow(M + 0.2e1 * r_H, 0.2e1) * (pow(M, 0.4e1) + 0.4607142857e1 * pow(M, 0.3e1) * r_H - 0.1205357143e2 * M * M * r_H * r_H - 0.1426190476e2 * pow(r_H, 0.3e1) * M + 0.4761904762e-1 * pow(r_H, 0.4e1)) * pow(x, 0.6e1) + 0.129e3 * pow(M + 0.2e1 * r_H, 0.2e1) * (M - 0.2e1 * r_H) * M * (M * M + 0.5581395349e1 * M * r_H + 0.3441860465e1 * r_H * r_H) * pow(x, 0.5e1) - 0.141e3 * pow(M + 0.2e1 * r_H, 0.2e1) * M * (pow(M, 0.3e1) + 0.2978723404e1 * M * M * r_H - 0.4446808511e1 * M * r_H * r_H - 0.3276595745e1 * pow(r_H, 0.3e1)) * pow(x, 0.4e1) + 0.114e3 * pow(M + 0.2e1 * r_H, 0.2e1) * M * (pow(M, 0.3e1) + 0.2526315789e1 * M * M * r_H - 0.2e1 * M * r_H * r_H - 0.1333333333e1 * pow(r_H, 0.3e1)) * pow(x, 0.3e1) - 0.66e2 * (pow(M, 0.3e1) + 0.2181818182e1 * M * M * r_H - 0.5454545455e0 * M * r_H * r_H - 0.3636363636e0 * pow(r_H, 0.3e1)) * pow(M + 0.2e1 * r_H, 0.2e1) * M * x * x + 0.24e2 * pow(M, 0.3e1) * pow(M + 0.2e1 * r_H, 0.3e1) * x - 0.4e1 * pow(M, 0.3e1) * pow(M + 0.2e1 * r_H, 0.3e1))) * r_H * r_H * M * pow(0.2500000000e0 * x * x * r_H * r_H * pow(x - 0.1e1, 0.2e1) * pow(x - 0.2e1, 0.2e1) * (M - 0.2e1 * r_H) * (M + 0.2e1 * r_H) * pow(cos(y), 0.2e1) + (-0.2500000000e0 * pow(x, 0.6e1) * pow(r_H, 0.3e1) + 0.7500000000e0 * r_H * r_H * (M + 0.2e1 * r_H) * pow(x, 0.5e1) - 0.1e1 * (M + 0.2e1 * r_H) * (M + 0.1750000000e1 * r_H) * r_H * pow(x, 0.4e1) + (M + 0.2e1 * r_H) * (M * M + 0.2e1 * M * r_H + 0.2e1 * r_H * r_H) * pow(x, 0.3e1) + (-0.3e1 * pow(M, 0.3e1) - 0.7e1 * M * M * r_H - 0.3e1 * M * r_H * r_H - 0.2e1 * pow(r_H, 0.3e1)) * x * x + 0.3e1 * M * M * (M + 0.2e1 * r_H) * x - 0.1e1 * M * M * (M + 0.2e1 * r_H)) * (-0.1e1 * x * x * r_H + (M + 0.2e1 * r_H) * x - 0.1e1 * M - 0.2e1 * r_H), -0.3e1));
}
#include <math.h>

double KERRd2u1dx2Xout (
  double x,
  double y,
  double r_H,
  double M)
{
  return(0.1562500000e-1 * (0.2e1 * pow(r_H, 0.4e1) * pow(x, 0.6e1) * pow(x - 0.1e1, 0.6e1) * pow(x - 0.2e1, 0.6e1) * pow(M - 0.2e1 * r_H, 0.4e1) * pow(M + 0.2e1 * r_H, 0.4e1) * pow(cos(y), 0.8e1) + 0.24e2 * pow(M + 0.2e1 * r_H, 0.3e1) * pow(x - 0.1e1, 0.4e1) * pow(M - 0.2e1 * r_H, 0.3e1) * r_H * r_H * (0.7500000000e0 * pow(x, 0.12e2) * pow(r_H, 0.4e1) - 0.1e1 * pow(r_H, 0.3e1) * (M + 0.9e1 * r_H) * pow(x, 0.11e2) + (0.4716666667e2 * pow(r_H, 0.4e1) + 0.1333333333e1 * M * M * r_H * r_H + 0.11e2 * pow(r_H, 0.3e1) * M) * pow(x, 0.10e2) + (-0.1416666667e3 * pow(r_H, 0.4e1) - 0.2e1 * pow(M, 0.3e1) * r_H - 0.1333333333e2 * M * M * r_H * r_H - 0.52e2 * pow(r_H, 0.3e1) * M) * pow(x, 0.9e1) + (pow(M, 0.4e1) + 0.18e2 * pow(M, 0.3e1) * r_H + 0.5866666667e2 * M * M * r_H * r_H + 0.138e3 * pow(r_H, 0.3e1) * M + 0.2683333333e3 * pow(r_H, 0.4e1)) * pow(x, 0.8e1) + (-0.8e1 * pow(M, 0.4e1) - 0.68e2 * pow(M, 0.3e1) * r_H - 0.1493333333e3 * M * M * r_H * r_H - 0.3306666667e3 * pow(r_H, 0.4e1) - 0.224e3 * pow(r_H, 0.3e1) * M) * pow(x, 0.7e1) + (0.24e2 * pow(M, 0.4e1) + 0.140e3 * pow(M, 0.3e1) * r_H + 0.2613333333e3 * pow(r_H, 0.4e1) + 0.2346666667e3 * M * M * r_H * r_H + 0.224e3 * pow(r_H, 0.3e1) * M) * pow(x, 0.6e1) + (-0.32e2 * pow(M, 0.4e1) - 0.1226666667e3 * pow(r_H, 0.4e1) - 0.154e3 * pow(M, 0.3e1) * r_H - 0.2133333333e3 * M * M * r_H * r_H - 0.128e3 * pow(r_H, 0.3e1) * M) * pow(x, 0.5e1) + (0.42e2 * pow(M, 0.3e1) * r_H + 0.5733333333e2 * M * M * r_H * r_H + 0.32e2 * pow(r_H, 0.3e1) * M + 0.2666666667e2 * pow(r_H, 0.4e1) + 0.9e1 * pow(M, 0.4e1)) * pow(x, 0.4e1) + 0.28e2 * M * M * pow(M + 0.2e1 * r_H, 0.2e1) * pow(x, 0.3e1) - 0.38e2 * M * M * pow(M + 0.2e1 * r_H, 0.2e1) * x * x + 0.20e2 * M * M * pow(M + 0.2e1 * r_H, 0.2e1) * x - 0.4e1 * M * M * pow(M + 0.2e1 * r_H, 0.2e1)) * pow(x - 0.2e1, 0.2e1) * x * x * pow(cos(y), 0.6e1) + 0.192e3 * pow(M + 0.2e1 * r_H, 0.2e1) * pow(x - 0.1e1, 0.2e1) * pow(M - 0.2e1 * r_H, 0.2e1) * (0.2187500000e0 * pow(x, 0.20e2) * pow(r_H, 0.8e1) - 0.1e1 * (0.4375000000e1 * r_H + M) * pow(r_H, 0.7e1) * pow(x, 0.19e2) + (0.4081250000e2 * pow(r_H, 0.8e1) + 0.2343750000e1 * M * M * pow(r_H, 0.6e1) + 0.19e2 * pow(r_H, 0.7e1) * M) * pow(x, 0.18e2) + (-0.2358750000e3 * pow(r_H, 0.8e1) - 0.3875000000e1 * pow(M, 0.3e1) * pow(r_H, 0.5e1) - 0.4218750000e2 * M * M * pow(r_H, 0.6e1) - 0.1667500000e3 * pow(r_H, 0.7e1) * M) * pow(x, 0.17e2) + (0.9455000000e3 * pow(r_H, 0.8e1) + 0.4145833333e1 * pow(M, 0.4e1) * pow(r_H, 0.4e1) + 0.6587500000e2 * pow(M, 0.3e1) * pow(r_H, 0.5e1) + 0.3506145833e3 * M * M * pow(r_H, 0.6e1) + 0.8967500000e3 * pow(r_H, 0.7e1) * M) * pow(x, 0.16e2) + (-0.2786e4 * pow(r_H, 0.8e1) - 0.3250000000e1 * pow(M, 0.5e1) * pow(r_H, 0.3e1) - 0.6633333333e2 * pow(M, 0.4e1) * pow(r_H, 0.4e1) - 0.5156250000e3 * pow(M, 0.3e1) * pow(r_H, 0.5e1) - 0.1784833333e4 * M * M * pow(r_H, 0.6e1) - 0.3303500000e4 * pow(r_H, 0.7e1) * M) * pow(x, 0.15e2) + 0.3125000000e1 * r_H * r_H * (pow(M, 0.6e1) + 0.1560000000e2 * pow(M, 0.5e1) * r_H + 0.1549200000e3 * r_H * r_H * pow(M, 0.4e1) + 0.7886000000e3 * pow(M, 0.3e1) * pow(r_H, 0.3e1) + 0.1990980000e4 * M * M * pow(r_H, 0.4e1) + 0.2822560000e4 * pow(r_H, 0.5e1) * M + 0.1994800000e4 * pow(r_H, 0.6e1)) * pow(x, 0.14e2) + (-0.1077650000e5 * pow(r_H, 0.8e1) - 0.2666666667e1 * pow(M, 0.7e1) * r_H - 0.4375000000e2 * pow(M, 0.6e1) * r_H * r_H - 0.3374166667e3 * pow(M, 0.5e1) * pow(r_H, 0.3e1) - 0.2134416667e4 * pow(M, 0.4e1) * pow(r_H, 0.4e1) - 0.8023e4 * pow(M, 0.3e1) * pow(r_H, 0.5e1) - 0.1571704167e5 * M * M * pow(r_H, 0.6e1) - 0.17595e5 * pow(r_H, 0.7e1) * M) * pow(x, 0.13e2) + (pow(M, 0.8e1) + 0.1448450000e5 * pow(r_H, 0.8e1) + 0.2916250000e3 * pow(M, 0.6e1) * r_H * r_H + 0.6327979167e4 * pow(M, 0.4e1) * pow(r_H, 0.4e1) + 0.2965379167e5 * M * M * pow(r_H, 0.6e1) + 0.3466666667e2 * pow(M, 0.7e1) * r_H + 0.1428916667e4 * pow(M, 0.5e1) * pow(r_H, 0.3e1) + 0.1873625000e5 * pow(M, 0.3e1) * pow(r_H, 0.5e1) + 0.26598e5 * pow(r_H, 0.7e1) * M) * pow(x, 0.12e2) + (-0.12e2 * pow(M, 0.8e1) - 0.15084e5 * pow(r_H, 0.8e1) - 0.1224500000e4 * pow(M, 0.6e1) * r_H * r_H - 0.1323675000e5 * pow(M, 0.4e1) * pow(r_H, 0.4e1) - 0.42318e5 * M * M * pow(r_H, 0.6e1) - 0.213e3 * pow(M, 0.7e1) * r_H - 0.4156750000e4 * pow(M, 0.5e1) * pow(r_H, 0.3e1) - 0.32032e5 * pow(M, 0.3e1) * pow(r_H, 0.5e1) - 0.30524e5 * pow(r_H, 0.7e1) * M) * pow(x, 0.11e2) + (0.6566666667e2 * pow(M, 0.8e1) + 0.11998e5 * pow(r_H, 0.8e1) + 0.3621625000e4 * pow(M, 0.6e1) * r_H * r_H + 0.1981758333e5 * pow(M, 0.4e1) * pow(r_H, 0.4e1) + 0.4543416667e5 * M * M * pow(r_H, 0.6e1) + 0.8176666667e3 * pow(M, 0.7e1) * r_H + 0.8877916667e4 * pow(M, 0.5e1) * pow(r_H, 0.3e1) + 0.39754e5 * pow(M, 0.3e1) * pow(r_H, 0.5e1) + 0.26312e5 * pow(r_H, 0.7e1) * M) * pow(x, 0.10e2) + (-0.2166666667e3 * pow(M, 0.8e1) - 0.7080e4 * pow(r_H, 0.8e1) - 0.8001250000e4 * pow(M, 0.6e1) * r_H * r_H - 0.2085033333e5 * pow(M, 0.4e1) * pow(r_H, 0.4e1) - 0.3538466667e5 * M * M * pow(r_H, 0.6e1) - 0.2179666667e4 * pow(M, 0.7e1) * r_H - 0.1461991667e5 * pow(M, 0.5e1) * pow(r_H, 0.3e1) - 0.33719e5 * pow(M, 0.3e1) * pow(r_H, 0.5e1) - 0.16592e5 * pow(r_H, 0.7e1) * M) * pow(x, 0.9e1) + (0.480e3 * pow(M, 0.8e1) + 0.2936e4 * pow(r_H, 0.8e1) + 0.1363987500e5 * pow(M, 0.6e1) * r_H * r_H + 0.14258e5 * pow(M, 0.4e1) * pow(r_H, 0.4e1) + 0.17516e5 * M * M * pow(r_H, 0.6e1) + 0.4239e4 * pow(M, 0.7e1) * r_H + 0.1939575000e5 * pow(M, 0.5e1) * pow(r_H, 0.3e1) + 0.14897e5 * pow(M, 0.3e1) * pow(r_H, 0.5e1) + 0.7264e4 * pow(r_H, 0.7e1) * M) * pow(x, 0.8e1) - 0.752e3 * (pow(M, 0.6e1) + 0.4183067376e1 * pow(M, 0.5e1) * r_H + 0.3471187943e1 * r_H * r_H * pow(M, 0.4e1) - 0.1967198582e1 * pow(M, 0.3e1) * pow(r_H, 0.3e1) + 0.1347517730e0 * M * M * pow(r_H, 0.4e1) + 0.4042553191e0 * pow(r_H, 0.5e1) * M + 0.2553191489e0 * pow(r_H, 0.6e1)) * pow(M + 0.2e1 * r_H, 0.2e1) * pow(x, 0.7e1) + 0.854e3 * pow(M + 0.2e1 * r_H, 0.2e1) * (pow(M, 0.6e1) + 0.3863778298e1 * pow(M, 0.5e1) * r_H + 0.2794594067e1 * r_H * r_H * pow(M, 0.4e1) - 0.2724043716e1 * pow(M, 0.3e1) * pow(r_H, 0.3e1) - 0.1662373146e1 * M * M * pow(r_H, 0.4e1) + 0.4683840749e-1 * pow(r_H, 0.5e1) * M + 0.2810304450e-1 * pow(r_H, 0.6e1)) * pow(x, 0.6e1) - 0.708e3 * pow(M + 0.2e1 * r_H, 0.3e1) * M * M * (pow(M, 0.3e1) + 0.1744350282e1 * M * M * r_H - 0.8425141243e0 * M * r_H * r_H - 0.9590395480e0 * pow(r_H, 0.3e1)) * pow(x, 0.5e1) + 0.425e3 * (pow(M, 0.3e1) + 0.1748235294e1 * M * M * r_H - 0.6717647059e0 * M * r_H * r_H - 0.8023529412e0 * pow(r_H, 0.3e1)) * pow(M + 0.2e1 * r_H, 0.3e1) * M * M * pow(x, 0.4e1) - 0.180e3 * pow(M + 0.2e1 * r_H, 0.3e1) * (pow(M, 0.3e1) + 0.1822222222e1 * M * M * r_H - 0.4444444444e0 * M * r_H * r_H - 0.5333333333e0 * pow(r_H, 0.3e1)) * M * M * pow(x, 0.3e1) + 0.51e2 * pow(M + 0.2e1 * r_H, 0.3e1) * M * M * (pow(M, 0.3e1) + 0.1921568627e1 * M * M * r_H - 0.1960784314e0 * M * r_H * r_H - 0.2352941176e0 * pow(r_H, 0.3e1)) * x * x - 0.8666666667e1 * pow(M, 0.4e1) * pow(M + 0.2e1 * r_H, 0.4e1) * x + 0.6666666667e0 * pow(M, 0.4e1) * pow(M + 0.2e1 * r_H, 0.4e1)) * pow(cos(y), 0.4e1) + 0.384e3 * (M + 0.2e1 * r_H) * (-0.9895833333e-1 * pow(x, 0.18e2) * pow(r_H, 0.7e1) + 0.5156250000e0 * (M + 0.3454545455e1 * r_H) * pow(r_H, 0.6e1) * pow(x, 0.17e2) + (-0.8765625000e1 * M * pow(r_H, 0.6e1) - 0.1312500000e1 * M * M * pow(r_H, 0.5e1) - 0.1468750000e2 * pow(r_H, 0.7e1)) * pow(x, 0.16e2) + (0.21e2 * M * M * pow(r_H, 0.5e1) + 0.2416666667e1 * pow(M, 0.3e1) * pow(r_H, 0.4e1) + 0.6728125000e2 * M * pow(r_H, 0.6e1) + 0.7350000000e2 * pow(r_H, 0.7e1)) * pow(x, 0.15e2) - 0.3062500000e1 * (pow(M, 0.4e1) + 0.1183673469e2 * pow(M, 0.3e1) * r_H + 0.4977551020e2 * M * M * r_H * r_H + 0.1005612245e3 * pow(r_H, 0.3e1) * M + 0.8134693878e2 * pow(r_H, 0.4e1)) * pow(r_H, 0.3e1) * pow(x, 0.14e2) + 0.3437500000e1 * (pow(M, 0.5e1) + 0.1247272727e2 * pow(M, 0.4e1) * r_H + 0.7152727273e2 * pow(M, 0.3e1) * r_H * r_H + 0.1932000000e3 * M * M * pow(r_H, 0.3e1) + 0.2719636364e3 * M * pow(r_H, 0.4e1) + 0.1756363636e3 * pow(r_H, 0.5e1)) * r_H * r_H * pow(x, 0.13e2) - 0.3e1 * (pow(M, 0.6e1) + 0.1489583333e2 * pow(M, 0.5e1) * r_H + 0.9080555556e2 * r_H * r_H * pow(M, 0.4e1) + 0.3324027778e3 * pow(M, 0.3e1) * pow(r_H, 0.3e1) + 0.6448750000e3 * M * M * pow(r_H, 0.4e1) + 0.6594791667e3 * pow(r_H, 0.5e1) * M + 0.3582638889e3 * pow(r_H, 0.6e1)) * r_H * pow(x, 0.12e2) + (M + 0.2e1 * r_H) * (pow(M, 0.6e1) + 0.34e2 * pow(M, 0.5e1) * r_H + 0.2025000000e3 * r_H * r_H * pow(M, 0.4e1) + 0.6345000000e3 * pow(M, 0.3e1) * pow(r_H, 0.3e1) + 0.1416125000e4 * M * M * pow(r_H, 0.4e1) + 0.1136750000e4 * pow(r_H, 0.5e1) * M + 0.7092500000e3 * pow(r_H, 0.6e1)) * pow(x, 0.11e2) - 0.11e2 * (pow(M, 0.6e1) + 0.1618181818e2 * pow(M, 0.5e1) * r_H + 0.5938636364e2 * r_H * r_H * pow(M, 0.4e1) + 0.1230738636e3 * pow(M, 0.3e1) * pow(r_H, 0.3e1) + 0.2094772727e3 * M * M * pow(r_H, 0.4e1) + 0.1145227273e3 * pow(r_H, 0.5e1) * M + 0.6270454545e2 * pow(r_H, 0.6e1)) * (M + 0.2e1 * r_H) * pow(x, 0.10e2) + 0.5466666667e2 * (M + 0.2e1 * r_H) * (pow(M, 0.6e1) + 0.1043902439e2 * pow(M, 0.5e1) * r_H + 0.2655373476e2 * r_H * r_H * pow(M, 0.4e1) + 0.3538871951e2 * pow(M, 0.3e1) * pow(r_H, 0.3e1) + 0.4826371951e2 * M * M * pow(r_H, 0.4e1) + 0.1771951220e2 * pow(r_H, 0.5e1) * M + 0.8829268293e1 * pow(r_H, 0.6e1)) * pow(x, 0.9e1) + (-0.162e3 * pow(M, 0.7e1) - 0.1574e4 * pow(M, 0.6e1) * r_H - 0.4856437500e4 * pow(M, 0.5e1) * r_H * r_H - 0.6485375000e4 * pow(M, 0.4e1) * pow(r_H, 0.3e1) - 0.5595e4 * pow(M, 0.3e1) * pow(r_H, 0.4e1) - 0.4595e4 * M * M * pow(r_H, 0.5e1) - 0.1221e4 * M * pow(r_H, 0.6e1) - 0.462e3 * pow(r_H, 0.7e1)) * pow(x, 0.8e1) + (0.318e3 * pow(M, 0.7e1) + 0.2608e4 * pow(M, 0.6e1) * r_H + 0.6842625000e4 * pow(M, 0.5e1) * r_H * r_H + 0.6616e4 * pow(M, 0.4e1) * pow(r_H, 0.3e1) + 0.2547500000e4 * pow(M, 0.3e1) * pow(r_H, 0.4e1) + 0.1972e4 * M * M * pow(r_H, 0.5e1) + 0.372e3 * M * pow(r_H, 0.6e1) + 0.136e3 * pow(r_H, 0.7e1)) * pow(x, 0.7e1) + (-0.434e3 * pow(M, 0.7e1) - 0.3164e4 * pow(M, 0.6e1) * r_H - 0.7360875000e4 * pow(M, 0.5e1) * r_H * r_H - 0.5339166667e4 * pow(M, 0.4e1) * pow(r_H, 0.3e1) + 0.4188333333e3 * pow(M, 0.3e1) * pow(r_H, 0.4e1) + 0.22e2 * M * M * pow(r_H, 0.5e1) - 0.52e2 * M * pow(r_H, 0.6e1) - 0.1866666667e2 * pow(r_H, 0.7e1)) * pow(x, 0.6e1) + 0.420e3 * pow(M + 0.2e1 * r_H, 0.2e1) * (pow(M, 0.3e1) + 0.2723809524e1 * M * M * r_H - 0.4940476190e0 * M * r_H * r_H - 0.4119047619e0 * pow(r_H, 0.3e1)) * M * M * pow(x, 0.5e1) - 0.288e3 * pow(M + 0.2e1 * r_H, 0.2e1) * (pow(M, 0.3e1) + 0.2371527778e1 * M * M * r_H - 0.4973958333e0 * M * r_H * r_H - 0.4288194444e0 * pow(r_H, 0.3e1)) * M * M * pow(x, 0.4e1) + 0.137e3 * pow(M + 0.2e1 * r_H, 0.2e1) * (pow(M, 0.3e1) + 0.2160583942e1 * M * M * r_H - 0.3576642336e0 * M * r_H * r_H - 0.3065693431e0 * pow(r_H, 0.3e1)) * M * M * pow(x, 0.3e1) - 0.43e2 * pow(M + 0.2e1 * r_H, 0.2e1) * M * M * (pow(M, 0.3e1) + 0.2046511628e1 * M * M * r_H - 0.1627906977e0 * M * r_H * r_H - 0.1395348837e0 * pow(r_H, 0.3e1)) * x * x + 0.8e1 * pow(M, 0.4e1) * pow(M + 0.2e1 * r_H, 0.3e1) * x - 0.6666666667e0 * pow(M, 0.4e1) * pow(M + 0.2e1 * r_H, 0.3e1)) * (M - 0.2e1 * r_H) * pow(-0.1e1 * x * x * r_H + (M + 0.2e1 * r_H) * x - 0.1e1 * M - 0.2e1 * r_H, 0.3e1) * pow(cos(y), 0.2e1) + 0.192e3 * (0.6250000000e-1 * pow(r_H, 0.8e1) * pow(x, 0.18e2) - 0.5000000000e0 * pow(r_H, 0.7e1) * (0.2250000000e1 * r_H + M) * pow(x, 0.17e2) + (0.9333333333e1 * pow(r_H, 0.8e1) + 0.1927083333e1 * M * M * pow(r_H, 0.6e1) + 0.8500000000e1 * pow(r_H, 0.7e1) * M) * pow(x, 0.16e2) - 0.4916666667e1 * (M * M + 0.4271186441e1 * M * r_H + 0.4813559322e1 * r_H * r_H) * (M + 0.2e1 * r_H) * pow(r_H, 0.5e1) * pow(x, 0.15e2) + (0.1639166667e3 * pow(r_H, 0.8e1) + 0.8927083333e1 * pow(M, 0.4e1) * pow(r_H, 0.4e1) + 0.7375000000e2 * pow(M, 0.3e1) * pow(r_H, 0.5e1) + 0.2233125000e3 * M * M * pow(r_H, 0.6e1) + 0.305e3 * pow(r_H, 0.7e1) * M) * pow(x, 0.14e2) - 0.1158333333e2 * (M + 0.2e1 * r_H) * pow(r_H, 0.3e1) * (pow(M, 0.4e1) + 0.8789568345e1 * pow(M, 0.3e1) * r_H + 0.2573021583e2 * M * M * r_H * r_H + 0.3211151079e2 * pow(r_H, 0.3e1) * M + 0.1767625899e2 * pow(r_H, 0.4e1)) * pow(x, 0.13e2) + 0.1054166667e2 * (M + 0.2e1 * r_H) * (pow(M, 0.5e1) + 0.1228458498e2 * pow(M, 0.4e1) * r_H + 0.5092193676e2 * pow(M, 0.3e1) * r_H * r_H + 0.9238537549e2 * M * M * pow(r_H, 0.3e1) + 0.8066403162e2 * M * pow(r_H, 0.4e1) + 0.3598418972e2 * pow(r_H, 0.5e1)) * r_H * r_H * pow(x, 0.12e2) + (-0.1052666667e4 * pow(r_H, 0.8e1) - 0.1265000000e3 * pow(M, 0.6e1) * r_H * r_H - 0.3050708333e4 * pow(M, 0.4e1) * pow(r_H, 0.4e1) - 0.5686e4 * M * M * pow(r_H, 0.6e1) - 0.5333333333e1 * pow(M, 0.7e1) * r_H - 0.9001666667e3 * pow(M, 0.5e1) * pow(r_H, 0.3e1) - 0.5580e4 * pow(M, 0.3e1) * pow(r_H, 0.5e1) - 0.3288e4 * pow(r_H, 0.7e1) * M) * pow(x, 0.11e2) + (M + 0.2e1 * r_H) * (pow(M, 0.7e1) + 0.5666666667e2 * pow(M, 0.6e1) * r_H + 0.5865000000e3 * pow(M, 0.5e1) * r_H * r_H + 0.2103166667e4 * pow(M, 0.4e1) * pow(r_H, 0.3e1) + 0.3648125000e4 * pow(M, 0.3e1) * pow(r_H, 0.4e1) + 0.3366416667e4 * M * M * pow(r_H, 0.5e1) + 0.1605333333e4 * M * pow(r_H, 0.6e1) + 0.544e3 * pow(r_H, 0.7e1)) * pow(x, 0.10e2) - 0.10e2 * pow(M + 0.2e1 * r_H, 0.2e1) * (pow(M, 0.6e1) + 0.2536666667e2 * pow(M, 0.5e1) * r_H + 0.1305333333e3 * r_H * r_H * pow(M, 0.4e1) + 0.1876916667e3 * pow(M, 0.3e1) * pow(r_H, 0.3e1) + 0.1615000000e3 * M * M * pow(r_H, 0.4e1) + 0.5553333333e2 * pow(r_H, 0.5e1) * M + 0.2046666667e2 * pow(r_H, 0.6e1)) * pow(x, 0.9e1) + 0.4466666667e2 * pow(M + 0.2e1 * r_H, 0.2e1) * (pow(M, 0.6e1) + 0.1576865672e2 * pow(M, 0.5e1) * r_H + 0.5388526119e2 * r_H * r_H * pow(M, 0.4e1) + 0.4568097015e2 * pow(M, 0.3e1) * pow(r_H, 0.3e1) + 0.3031343284e2 * M * M * pow(r_H, 0.4e1) + 0.6925373134e1 * pow(r_H, 0.5e1) * M + 0.2388059701e1 * pow(r_H, 0.6e1)) * pow(x, 0.8e1) - 0.1173333333e3 * pow(M + 0.2e1 * r_H, 0.2e1) * (pow(M, 0.6e1) + 0.1109659091e2 * pow(M, 0.5e1) * r_H + 0.2701420455e2 * r_H * r_H * pow(M, 0.4e1) + 0.1193181818e2 * pow(M, 0.3e1) * pow(r_H, 0.3e1) + 0.6244318182e1 * M * M * pow(r_H, 0.4e1) + 0.8863636364e0 * pow(r_H, 0.5e1) * M + 0.2954545455e0 * pow(r_H, 0.6e1)) * pow(x, 0.7e1) + 0.2006666667e3 * pow(M + 0.2e1 * r_H, 0.2e1) * (pow(M, 0.6e1) + 0.8395348837e1 * pow(M, 0.5e1) * r_H + 0.1537001661e2 * r_H * r_H * pow(M, 0.4e1) + 0.2237541528e1 * pow(M, 0.3e1) * pow(r_H, 0.3e1) + 0.7890365448e0 * M * M * pow(r_H, 0.4e1) + 0.7973421927e-1 * pow(r_H, 0.5e1) * M + 0.2657807309e-1 * pow(r_H, 0.6e1)) * pow(x, 0.6e1) - 0.2333333333e3 * pow(M + 0.2e1 * r_H, 0.2e1) * M * M * (pow(M, 0.4e1) + 0.6681428571e1 * pow(M, 0.3e1) * r_H + 0.9695000000e1 * M * M * r_H * r_H - 0.6000000000e0 * pow(r_H, 0.3e1) * M - 0.4714285714e0 * pow(r_H, 0.4e1)) * pow(x, 0.5e1) + 0.1866666667e3 * pow(M + 0.2e1 * r_H, 0.2e1) * (pow(M, 0.4e1) + 0.5537500000e1 * pow(M, 0.3e1) * r_H + 0.6743750000e1 * M * M * r_H * r_H - 0.1178571429e1 * pow(r_H, 0.3e1) * M - 0.6321428571e0 * pow(r_H, 0.4e1)) * M * M * pow(x, 0.4e1) - 0.1013333333e3 * pow(M + 0.2e1 * r_H, 0.3e1) * M * M * (pow(M, 0.3e1) + 0.2763157895e1 * M * M * r_H - 0.3552631579e0 * M * r_H * r_H - 0.2368421053e0 * pow(r_H, 0.3e1)) * pow(x, 0.3e1) + 0.3566666667e2 * pow(M + 0.2e1 * r_H, 0.3e1) * M * M * (pow(M, 0.3e1) + 0.2261682243e1 * M * M * r_H - 0.1682242991e0 * M * r_H * r_H - 0.1121495327e0 * pow(r_H, 0.3e1)) * x * x - 0.7333333333e1 * pow(M, 0.4e1) * pow(M + 0.2e1 * r_H, 0.4e1) * x + 0.6666666667e0 * pow(M, 0.4e1) * pow(M + 0.2e1 * r_H, 0.4e1)) * pow(-0.1e1 * x * x * r_H + (M + 0.2e1 * r_H) * x - 0.1e1 * M - 0.2e1 * r_H, 0.4e1)) * pow(0.2500000000e0 * x * x * r_H * r_H * pow(x - 0.1e1, 0.2e1) * pow(x - 0.2e1, 0.2e1) * (M - 0.2e1 * r_H) * (M + 0.2e1 * r_H) * pow(cos(y), 0.2e1) + (-0.2500000000e0 * pow(x, 0.6e1) * pow(r_H, 0.3e1) + 0.7500000000e0 * r_H * r_H * (M + 0.2e1 * r_H) * pow(x, 0.5e1) - 0.1e1 * (M + 0.2e1 * r_H) * (M + 0.1750000000e1 * r_H) * r_H * pow(x, 0.4e1) + (M + 0.2e1 * r_H) * (M * M + 0.2e1 * M * r_H + 0.2e1 * r_H * r_H) * pow(x, 0.3e1) + (-0.3e1 * pow(M, 0.3e1) - 0.7e1 * M * M * r_H - 0.3e1 * M * r_H * r_H - 0.2e1 * pow(r_H, 0.3e1)) * x * x + 0.3e1 * M * M * (M + 0.2e1 * r_H) * x - 0.1e1 * M * M * (M + 0.2e1 * r_H)) * (-0.1e1 * x * x * r_H + (M + 0.2e1 * r_H) * x - 0.1e1 * M - 0.2e1 * r_H), -0.3e1));
}
double KERRd2u2dx2Xout (
  double x,
  double y,
  double r_H,
  double M)
{
  return(0.12e2 * x * x - 0.24e2 * x + 0.8e1);
}
#include <math.h>

double KERRd2u3dx2Xout (
  double x,
  double y,
  double r_H,
  double M)
{
  return(-0.1e1 * (0.3750000000e0 * pow(M + 0.2e1 * r_H, 0.2e1) * pow(x - 0.1e1, 0.6e1) * pow(M - 0.2e1 * r_H, 0.2e1) * ((-0.5000000000e0 * pow(x, 0.4e1) + 0.2e1 * pow(x, 0.3e1) - 0.6e1 * x * x + 0.8e1 * x - 0.4e1) * r_H + M * (x - 0.1e1) * (x * x - 0.2e1 * x + 0.2e1)) * pow(r_H, 0.4e1) * pow(x - 0.2e1, 0.2e1) * x * x * pow(cos(y), 0.4e1) + 0.7500000000e0 * (M + 0.2e1 * r_H) * pow(x - 0.1e1, 0.2e1) * (M - 0.2e1 * r_H) * r_H * r_H * (-0.7500000000e0 * (pow(x, 0.8e1) - 0.8e1 * pow(x, 0.7e1) + 0.3333333333e2 * pow(x, 0.6e1) - 0.88e2 * pow(x, 0.5e1) + 0.152e3 * pow(x, 0.4e1) - 0.1706666667e3 * pow(x, 0.3e1) + 0.1226666667e3 * x * x - 0.5333333333e2 * x + 0.1066666667e2) * (x * x - 0.2e1 * x + 0.2e1) * pow(x - 0.2e1, 0.2e1) * x * x * pow(r_H, 0.5e1) + 0.4083333333e1 * (x - 0.1e1) * (pow(x, 0.8e1) - 0.8e1 * pow(x, 0.7e1) + 0.3220408163e2 * pow(x, 0.6e1) - 0.8122448980e2 * pow(x, 0.5e1) + 0.1362448980e3 * pow(x, 0.4e1) - 0.1528163265e3 * pow(x, 0.3e1) + 0.1116734694e3 * x * x - 0.4897959184e2 * x + 0.9795918367e1) * M * pow(x - 0.2e1, 0.2e1) * x * x * pow(r_H, 0.4e1) - 0.8500000000e1 * pow(x - 0.1e1, 0.4e1) * (pow(x, 0.8e1) - 0.8e1 * pow(x, 0.7e1) + 0.2901960784e2 * pow(x, 0.6e1) - 0.6211764706e2 * pow(x, 0.5e1) + 0.8235294118e2 * pow(x, 0.4e1) - 0.6462745098e2 * pow(x, 0.3e1) + 0.2447058824e2 * x * x - 0.1254901961e1) * M * M * pow(r_H, 0.3e1) + 0.9e1 * (pow(x, 0.6e1) - 0.6e1 * pow(x, 0.5e1) + 0.1533333333e2 * pow(x, 0.4e1) - 0.2133333333e2 * pow(x, 0.3e1) + 0.1422222222e2 * x * x - 0.1777777778e1 * x - 0.1777777778e1) * pow(x - 0.1e1, 0.5e1) * pow(M, 0.3e1) * r_H * r_H - 0.5e1 * (pow(x, 0.4e1) - 0.4e1 * pow(x, 0.3e1) + 0.4400000000e1 * x * x - 0.8000000000e0 * x - 0.8000000000e0) * pow(x - 0.1e1, 0.4e1) * (x * x - 0.2e1 * x + 0.2e1) * pow(M, 0.4e1) * r_H + pow(x - 0.1e1, 0.5e1) * pow(M, 0.5e1) * (pow(x, 0.4e1) - 0.4e1 * pow(x, 0.3e1) + 0.4666666667e1 * x * x - 0.1333333333e1 * x - 0.1333333333e1)) * pow(cos(y), 0.2e1) + (0.6250000000e0 * (pow(x, 0.8e1) - 0.8e1 * pow(x, 0.7e1) + 0.2910000000e2 * pow(x, 0.6e1) - 0.6260000000e2 * pow(x, 0.5e1) + 0.8860000000e2 * pow(x, 0.4e1) - 0.8640000000e2 * pow(x, 0.3e1) + 0.5760000000e2 * x * x - 0.24e2 * x + 0.4800000000e1) * pow(x - 0.2e1, 0.2e1) * x * x * pow(r_H, 0.6e1) - 0.2812500000e1 * pow(x - 0.1e1, 0.3e1) * (pow(x, 0.4e1) - 0.4e1 * pow(x, 0.3e1) + 0.7200000000e1 * x * x - 0.6400000000e1 * x + 0.3200000000e1) * M * pow(x - 0.2e1, 0.2e1) * x * x * pow(r_H, 0.5e1) + 0.5625000000e1 * pow(x - 0.1e1, 0.2e1) * (pow(x, 0.8e1) - 0.8e1 * pow(x, 0.7e1) + 0.2686666667e2 * pow(x, 0.6e1) - 0.4920000000e2 * pow(x, 0.5e1) + 0.5248888889e2 * pow(x, 0.4e1) - 0.3128888889e2 * pow(x, 0.3e1) + 0.8e1 * x * x + 0.7111111111e0 * x - 0.7111111111e0) * M * M * pow(r_H, 0.4e1) - 0.6500000000e1 * pow(x - 0.1e1, 0.3e1) * (pow(x, 0.6e1) - 0.6e1 * pow(x, 0.5e1) + 0.1384615385e2 * pow(x, 0.4e1) - 0.1538461538e2 * pow(x, 0.3e1) + 0.7384615385e1 * x * x - 0.1230769231e1) * pow(M, 0.3e1) * pow(r_H, 0.3e1) + 0.5250000000e1 * pow(x - 0.1e1, 0.4e1) * (pow(x, 0.4e1) - 0.4e1 * pow(x, 0.3e1) + 0.5238095238e1 * x * x - 0.2476190476e1 * x - 0.1904761905e0) * pow(M, 0.4e1) * r_H * r_H - 0.3e1 * pow(M, 0.5e1) * pow(x - 0.1e1, 0.7e1) * r_H + pow(M, 0.6e1) * pow(x - 0.1e1, 0.6e1)) * pow((-0.1e1 * x * x + 0.2e1 * x - 0.2e1) * r_H + M * (x - 0.1e1), 0.3e1)) * sqrt(M * M - 0.4e1 * r_H * r_H) * r_H * M * pow(0.2500000000e0 * x * x * r_H * r_H * pow(x - 0.1e1, 0.2e1) * pow(x - 0.2e1, 0.2e1) * (M - 0.2e1 * r_H) * (M + 0.2e1 * r_H) * pow(cos(y), 0.2e1) + (-0.2500000000e0 * x * x * (x * x - 0.2e1 * x + 0.2e1) * pow(x - 0.2e1, 0.2e1) * pow(r_H, 0.3e1) + 0.7500000000e0 * x * x * M * (x - 0.1e1) * pow(x - 0.2e1, 0.2e1) * r_H * r_H - 0.1e1 * M * M * (x * x - 0.2e1 * x + 0.2e1) * pow(x - 0.1e1, 0.2e1) * r_H + pow(x - 0.1e1, 0.3e1) * pow(M, 0.3e1)) * ((-0.1e1 * x * x + 0.2e1 * x - 0.2e1) * r_H + M * (x - 0.1e1)), -0.3e1));
}
double KERRd2u4dx2Xout (
  double x,
  double y,
  double r_H,
  double M)
{
  return(0.0e0);
}
