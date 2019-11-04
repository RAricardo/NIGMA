#include "include/nonlinear.h"
#include <ctype.h>
#include <string.h>
#include <jni.h>

extern "C"
JNIEXPORT void JNICALL
Java_com_example_nigma_JNICaller_setf(JNIEnv *env, jclass type, jstring expr_) {
    const char *expr = env->GetStringUTFChars(expr_, 0);

    setF(expr);

    env->ReleaseStringUTFChars(expr_, expr);
}

extern "C"
JNIEXPORT void JNICALL
Java_com_example_nigma_JNICaller_setg(JNIEnv *env, jclass type, jstring expr_) {
    const char *expr = env->GetStringUTFChars(expr_, 0);

    setG(expr);

    env->ReleaseStringUTFChars(expr_, expr);
}

extern "C"
JNIEXPORT void JNICALL
Java_com_example_nigma_JNICaller_setfdx(JNIEnv *env, jclass type, jstring expr_) {
    const char *expr = env->GetStringUTFChars(expr_, 0);

    setDervF(expr);

    env->ReleaseStringUTFChars(expr_, expr);
}

extern "C"
JNIEXPORT void JNICALL
Java_com_example_nigma_JNICaller_setfdx2(JNIEnv *env, jclass type, jstring expr_) {
    const char *expr = env->GetStringUTFChars(expr_, 0);

    setDervSecF(expr);

    env->ReleaseStringUTFChars(expr_, expr);
}

extern "C"
JNIEXPORT jstring JNICALL
Java_com_example_nigma_JNICaller_incrementalSearch(JNIEnv *env, jclass type, jdouble x0,
                                                   jdouble delta, jint niter) {
    return env->NewStringUTF(incrementalSearch(x0, delta, (int)niter).c_str());
}

extern "C"
JNIEXPORT jstring JNICALL
Java_com_example_nigma_JNICaller_bisection(JNIEnv *env, jclass type, jdouble xi, jdouble xs,
                                           jdouble tol, jint niter, jboolean error) {
    return env->NewStringUTF(bisection(xi, xs, tol, (int) niter, error).c_str());
}

extern "C"
JNIEXPORT jstring JNICALL
Java_com_example_nigma_JNICaller_regulaFalsi(JNIEnv *env, jclass type, jdouble xi, jdouble xs,
                                             jdouble tol, jint niter, jboolean error) {
    return env->NewStringUTF(regulaFalsi(xi, xs, tol, (int) niter, error).c_str());
}

extern "C"
JNIEXPORT jstring JNICALL
Java_com_example_nigma_JNICaller_fixedPoint(JNIEnv *env, jclass type, jdouble tol, jdouble xa,
                                            jdouble niter, jboolean relativeErr) {
    return env->NewStringUTF(fixedPoint(tol, xa, niter, relativeErr).c_str());
}

extern "C"
JNIEXPORT jstring JNICALL
Java_com_example_nigma_JNICaller_newtonMethod(JNIEnv *env, jclass type, jdouble tol, jdouble x0,
                                              jdouble niter, jboolean relativeErr) {
    return env->NewStringUTF(newtonMethod(tol,x0,niter,relativeErr).c_str());
}

extern "C"
JNIEXPORT jstring JNICALL
Java_com_example_nigma_JNICaller_secantMethod(JNIEnv *env, jclass type, jdouble tol, jdouble x0,
                                              jdouble x1, jdouble niter, jboolean relativeErr) {
    return env->NewStringUTF(secantMethod(tol,x0,x1,niter,relativeErr).c_str());
}

extern "C"
JNIEXPORT jstring JNICALL
Java_com_example_nigma_JNICaller_mpMethod(JNIEnv *env, jclass type, jdouble tol, jdouble x0,
                                          jdouble niter, jboolean relativeErr) {
    return env->NewStringUTF(multipleRootsMethod(tol,x0,niter,relativeErr).c_str());
}extern "C"
JNIEXPORT void JNICALL
Java_com_example_nigma_JNICaller_prepareCall(JNIEnv *env, jclass clazz) {
    prepareCall();
}