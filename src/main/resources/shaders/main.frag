#version 430
in vec3 Normal;
in vec3 FragPos;
in vec3 baseColor;
out vec4 color;
uniform mat4 modelMat;
uniform mat4 viewMat;
uniform mat4 mv_matrix;
uniform mat4 projMat;
uniform mat4 mvInverse;
uniform vec3 col;
uniform float alpha;
uniform vec3 lightColor;
uniform vec3 lightPos;
uniform float ambientStrength;
void main(void)
{
    vec3 ambient = ambientStrength * baseColor;
    vec3 normal2 = (mvInverse * vec4(Normal,0.0)).xyz;

    vec3 lightDir = normalize(lightPos - FragPos);
    float diff = max(dot(lightDir, normal2), 0.0f);
    vec3 diffuse = baseColor * lightColor * diff;
    //vec3 result = (ambient + diffuse) * baseColor;
    vec3 result = ambient + diffuse;
    color = vec4(result, alpha);
}