/*

 Author: Yannis STEFANELLI

 Creation Date: 01-06-2023 15:10:28

 Description :

*/

uniform float iTime;
uniform vec2 iResolution;
uniform vec2 iMouse;

vec2 rotate(vec2 uv, float th) {
  return mat2(cos(th), sin(th), -sin(th), cos(th)) * uv;
}

vec3 circle(vec2 uv, float r) {

    float y = length(vec2(uv.x, uv.y)) - r;

    return y > 0. ? vec3(0.) : 0.5 + 0.5 * cos(iTime + uv.xyx + vec3(0,2,4));
  
}

vec3 square(vec2 uv, float size, vec2 offset) {
    float x = uv.x - offset.x;
    float y = uv.y - offset.y;

    vec2 turn = rotate(vec2(x,y), iTime + sin(iTime));
    float d = max(abs(turn.x), abs(turn.y)) - size;
  
    return d > 0. ? vec3(1.) : 0.5 + 0.5 * cos(iTime + uv.xyx + vec3(1., 0., 0.));
}

void main()
{
    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = .5 * gl_FragCoord/iResolution.xy;
    uv -= .5;
    uv.x *= iResolution.x/iResolution.y;

    // Time varying pixel color
    vec3 col = square(uv, .2, 0);

    // col = vec3(uv.x, uv.y, 0);

    // Output to screen
    gl_FragColor = vec4(col,1.0);
}