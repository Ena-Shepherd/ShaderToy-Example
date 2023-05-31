
uniform float iTime;
uniform vec2 iResolution;

mat2 rotate(float a) {
  float s = sin(a);
  float c = cos(a);
  return mat2(c, -s, s, c);
}

void main()
{
    
    float t = iTime;
    vec4 o, ot, oc;

    vec3 p = vec3(gl_FragCoord.xy*1.-iResolution.xy,0);
    p.xy *= rotate(t*.2);
    
    p = normalize(p)*pow(length(p*.1),.2);
    p.z -= t*.5;
    float c = 0.;
    for(float a=.5; a>.1; a/=2.) 
        c += abs(dot(sin(p/a), cos(p.yzx/a)))*a;
        
    oc = vec4(abs(c-p.x*.5),abs(c-p.y*.5),c,1.);
    ot = vec4(step(c,.9));
    if (ot.r == 1.) {
        o = oc;
        o *= .1;
    }
    ot = vec4(step(c,.7));
    if (ot.r == 1.) {
        o = oc;
        o *= .3;
    }
    ot = vec4(step(c,.5));
    if (ot.r == 1.) {
        o = oc;
        o *= .9;
    }

    // Output to screen
    gl_FragColor = o;
}