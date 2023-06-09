/*

 Author: Yannis STEFANELLI

 Creation Date: 30-05-2023 13:38:46

 Description :

*/

uniform float iTime;
uniform vec2 iResolution;
uniform vec2 iMouse;

float distLine(vec3 ro, vec3 rd, vec3 pnt) {
    
    return length(cross(pnt - ro, rd)) / length(rd);
}

void main()
{
    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = gl_FragCoord/iResolution.xy;
    uv -= .5;
    
    uv.x *= iResolution.x / iResolution.y;
    
    
    // Time varying pixel color
    vec3 ro = vec3(0., 0., -2.);
    vec3 rd = vec3(uv, 0.) - ro;
    
    float t = iTime;
    vec3 p = vec3(sin(t), 0., 0.3+cos(t)); 
    float d = distLine(ro, rd, p);
    
    d = smoothstep(.1, .002, d);
    // Output to screen
    gl_FragColor = vec4(d);
}