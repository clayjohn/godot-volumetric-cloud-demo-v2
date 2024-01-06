shader_type spatial;

uniform sampler2D heightmap : filter_linear_mipmap;
uniform vec2 res = vec2(3416.0, 3834.0);
uniform vec4 ground : source_color;

varying vec2 coord;
varying float height;

//https://atyuwen.github.io/posts/normal-reconstruction
// https://www.shadertoy.com/view/fsVczR

vec3 getPos(vec2 fragCoord, float depth )
{
    vec2 p = fragCoord * res;
	return vec3(p.x, depth * 500.0, p.y);
}
vec3 computeNormalImproved(sampler2D depth, vec2 p, float mip)
{
    float c0 = textureLod(heightmap,p            , mip).r;
    float l1 = textureLod(depth,p-vec2(1,0) / res, mip).r;
    float l2 = textureLod(depth,p-vec2(2,0) / res, mip).r;
    float r1 = textureLod(depth,p+vec2(1,0) / res, mip).r;
    float r2 = textureLod(depth,p+vec2(2,0) / res, mip).r;
    float b2 = textureLod(depth,p-vec2(0,2) / res, mip).r;
    float b1 = textureLod(depth,p-vec2(0,1) / res, mip).r;
    float t1 = textureLod(depth,p+vec2(0,1) / res, mip).r;
    float t2 = textureLod(depth,p+vec2(0,2) / res, mip).r;
    
    float dl = abs(l1*l2/(2.0*l2-l1)-c0);
    float dr = abs(r1*r2/(2.0*r2-r1)-c0);
    float db = abs(b1*b2/(2.0*b2-b1)-c0);
    float dt = abs(t1*t2/(2.0*t2-t1)-c0);
    
    vec3 ce = getPos(p,c0);

    vec3 dpdx = (dl<dr) ?  ce-getPos(p-vec2(1,0) / res,l1) : 
                          -ce+getPos(p+vec2(1,0) / res,r1) ;
    vec3 dpdy = (db<dt) ?  ce-getPos(p-vec2(0,1) / res,b1) : 
                          -ce+getPos(p+vec2(0,1) / res,t1) ;

    return normalize(cross(dpdx,dpdy));
	//return vec3(c0);
}

void vertex() {
	coord = VERTEX.xz * 0.1 + 0.5;
	height = texture(heightmap, coord).r;
	VERTEX.y += height;
	VERTEX.x *= res.x / res.y;
}

void fragment() {
	ALBEDO = ground.rgb;
	vec3 norm = computeNormalImproved(heightmap, coord, 0.0);
	NORMAL = mat3(VIEW_MATRIX) * -norm;
}
