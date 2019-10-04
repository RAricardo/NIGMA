#include "include/nonlinear.h"
#include <ctype.h>
#include <string.h>
#include <jni.h>

/*JNIEXPORT jdouble JNICALL
Java_com_example_nigma2_MainActivity_f(JNIEnv *env, jobject instance, jdouble x) {
    return f(x);
}

JNIEXPORT void JNICALL
Java_com_example_nigma2_MainActivity_setf(JNIEnv *env, jobject instance, jstring expr_) {
    const char *expr = env->GetStringUTFChars(expr_, 0);

    setF(expr);

    env->ReleaseStringUTFChars(expr_, expr);
}

JNIEXPORT jstring JNICALL
Java_com_example_nigma2_MainActivity_incrementalSearch
(JNIEnv *env, jobject instance, jdouble x0, jdouble delta, jint niter){

}*/
/*
extern "C"
JNIEXPORT void JNICALL
Java_com_example_nigma_MainActivity_setf(JNIEnv *env, jobject instance, jstring expr_) {
    const char *expr = env->GetStringUTFChars(expr_, 0);

    setF(expr);

    env->ReleaseStringUTFChars(expr_, expr);
}

extern "C"
JNIEXPORT jstring JNICALL
Java_com_example_nigma_MainActivity_incrementalSearch(JNIEnv *env, jobject instance, jdouble x0,
                                                      jdouble delta, jint niter) {
    return env->NewStringUTF(incrementalSearch(x0, delta, (int)niter).c_str());}*/
extern "C"
JNIEXPORT void JNICALL
Java_com_example_nigma_IncrementalSearch_setf(JNIEnv *env, jobject instance, jstring expr_) {
    const char *expr = env->GetStringUTFChars(expr_, 0);

    setF(expr);

    env->ReleaseStringUTFChars(expr_, expr);
}extern "C"
JNIEXPORT jstring JNICALL
Java_com_example_nigma_IncrementalSearch_incrementalSearch(JNIEnv *env, jobject instance,
                                                           jdouble x0, jdouble delta, jint niter) {
    return env->NewStringUTF(incrementalSearch(x0, delta, (int)niter).c_str());
}