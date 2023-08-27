# OpenGLStarter

[Back to HOME](../index.md)

## C++ Easing Equations

Look at the example below:

```cpp
//with numbers
{
  float a = 0.0f;
  float b = 1.0f;
  float c = easeOutExpo(a,b,0.5f);
}
//with vectors
{
  vec3 a(1.0f,0.0f,0.0f);
  vec3 b(0.0f,1.0f,0.0f);
  vec3 c = easeOutElastic(a,b,0.5f);
}
//with mat4
{
  mat4 a = translate(vec3(1.0f,0.0f,0.0f));
  mat4 b = translate(vec3(0.0f, 1.0f, 0.0f));
  mat4 c = easeOutCirc(a,b,0.5f);
}
```

The project __aribeiro-text-easing-eq__ render all equations. You can see the print screen below:

![Robert Penner's equations](../images/easing-eq.png)
