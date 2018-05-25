#include <iostream>
#include "glm/glm.hpp"
#include "SDL.h"
#include "SDLauxiliary.h"
#include "TestModel.h"
#include <math.h> //For pc stuff
#include <algorithm> //pc stuff


#define M_PI 3.14159265358979323846f

using namespace std;
using glm::vec3;
using glm::ivec2;
using glm::mat3;

using glm::vec2;

// ----------------------------------------------------------------------------
// GLOBAL VARIABLES
float focalLength = 500; //Just field of view, see the entire image if the same as W and H
const int SCREEN_WIDTH = 500;
const int SCREEN_HEIGHT = 500;
SDL_Surface* screen;
int t;
vector<Triangle> triangles;
vec3 cameraPos(0, 0, -3.001); // (moved in cameraX, moved in cameraY, distance)!WRONG! What the camera sees as the images origin
mat3 R;
mat3 Ry;//rotation matrix y
mat3 Rx;//rotation matrix x
float yaw = 0.0f; // Yaw angle controlling camera rotation around y-axis
float xaw = 0.0f; // xaw angle controlling camera rotation around x - axis
vec3 currentColor;
//5
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH]; //inverse depth 1=z

												//6
vec3 lightPos(0, -0.5, -0.7);
vec3 lightPower = 14.1f * vec3(1, 1, 1);  //Increased lightpower because of bad lighting
vec3 indirectLightPowerPerArea = 0.5f*vec3(1, 1, 1);
vec3 currentNormal;
vec3 currentReflectance;

// ----------------------------------------------------------------------------
// STRUCTS

struct Pixel {
	int x;
	int y;
	float zinv;
	//vec3 illumination; //Not used anymore
	vec3 pos3d;
};
//6

struct Vertex
{
	vec3 position;
	//vec3 normal; //Not used anymore
	//vec3 reflectance; //Not used anymore
};
// ----------------------------------------------------------------------------
// FUNCTIONS

void Update();
void Draw();
/* OLD vertexshader
void VertexShader(const Vertex& v, Pixel& p) { //Fungerar som det ska
{
vec3 pos = (v.position - cameraPos )*R;
p.zinv = 1.0f / pos.z;
p.x = focalLength * pos.x * p.zinv + SCREEN_WIDTH / 2;
p.y = focalLength * pos.y * p.zinv + SCREEN_HEIGHT / 2;
}

//step 6 implementing eq 10 & 11
float r = glm::length(v.position - cameraPos);
vec3 specialr = glm::normalize(lightPos - v.position); //vector for shadowbeam (?) direction from surface -> light //p.pos3d = v.position
vec3 normal = v.normal;
float inner = glm::dot(normal, specialr);
vec3 D = lightPower * SDL_max(inner,0) / (4  * M_PI *r*r);
p.illumination = D * v.reflectance + indirectLightPowerPerArea;
}
*/

void VertexShader(const Vertex& v, Pixel& p) { //Instructions
	{
		vec3 pos = (v.position - cameraPos)*R;
		p.zinv = 1.0f / pos.z; //1/Z
		p.x = focalLength * pos.x * p.zinv + SCREEN_WIDTH / 2;
		p.y = focalLength * pos.y * p.zinv + SCREEN_HEIGHT / 2;
		p.pos3d = v.position; //give pixel location
	}

	/*//step 6 implementing eq 10 & 11
	float r = glm::length(v.position - lightPos); //RAAAAAAAAAAAAAAAAA! v.position - cameraPos e fel, sug min kuk
	vec3 specialr = glm::normalize(lightPos - v.position); //vector for shadowbeam (?) direction from surface -> light //p.pos3d = v.position
	vec3 normal = v.normal;
	float inner = glm::dot(normal, specialr);
	vec3 D = lightPower * SDL_max(inner, 0) / (4*M_PI*r*r);
	p.illumination = v.reflectance * (D + indirectLightPowerPerArea);
	*/
}




void PixelShader(const Pixel& p) //från Pek
{
	int x = p.x;
	int y = p.y;
	if (p.zinv > depthBuffer[y][x])
	{  // Basically the same as in the old vertexshader,implementing eq 10 & 11
		depthBuffer[y][x] = p.zinv; //corrected f. to p.
		float r = glm::length(p.pos3d - lightPos);
		vec3 specialr = glm::normalize(lightPos - p.pos3d); //vector for shadowray (?) direction from surface -> light //p.pos3d = v.position
		vec3 normal = currentNormal;
		float inner = glm::dot(normal, specialr);
		vec3 D = lightPower * SDL_max(inner, 0) / (4 * M_PI*r*r);
		vec3 currentLight = currentReflectance * (D + indirectLightPowerPerArea);
		PutPixelSDL(screen, x, y, currentLight); /* * currentColor*/

	}
}


void Interpolate(Pixel a, Pixel b, vector<Pixel>& result) //Fix comments and stuff, ändra runt
{
	int N = result.size(); //value a and the end value b and fills the std::vector result with values linearly interpolated in between.
	float maxN = (glm::max(N - 1, 1));
	float StepinX = (b.x - a.x) / maxN;
	float StepinY = (b.y - a.y) / maxN;
	float StepinZ = (b.zinv - a.zinv) / maxN;

	//Interpolation convertion, interpolate the transformed position (which is multiplied by inverse z value)

	b.pos3d = (b.pos3d * b.zinv); // b * 1/z according to hints
	a.pos3d = (a.pos3d * a.zinv); // a * 1/z

	vec3 Step3D = (b.pos3d - a.pos3d) / maxN; //adding for 6.1 //6.2 change all positions to pos3d
	Pixel current(a);
	for (int i = 0; i<N; ++i)
	{
		result[i].x = a.x + i * StepinX;
		result[i].y = a.y + i * StepinY;
		result[i].zinv = a.zinv + i * StepinZ;
		result[i].pos3d = a.pos3d + (float)i * Step3D;
		result[i].pos3d = (result[i].pos3d / result[i].zinv); //Restore by being divided by the zinverse

	}
}



void DrawLineSDL(SDL_Surface* surface, Pixel a, Pixel b, vec3 color) //Lab p.5
{
	Pixel delta;
	delta.x = glm::abs(a.x - b.x);
	delta.y = glm::abs(a.y - b.y);

	int pixels = glm::max(delta.x, delta.y) + 1;

	vector<Pixel> line(pixels);
	Interpolate(a, b, line);

	for (int i = 0; i < line.size(); ++i) 
	{
		//if (line[i].x >= 0 && line[i].x < SCREEN_WIDTH && line[i].y >= 0 && line[i].y < SCREEN_HEIGHT) //ha kvar??
		//{
			PixelShader(line[i]);
		//}
	}
}


/* Används inte längre
void DrawPolygonEdges(const vector<vec3>& vertices) {
int V = vertices.size();

// Transform each vertex from 3D world position to 2D image position:
vector<Pixel> projectedVertices(V);
for (int i = 0; i<V; ++i) {
VertexShader(vertices[i], projectedVertices[i]);
}

// Loop over all vertices and draw the edge from it to the next vertex:
for (int i = 0; i<V; ++i) {
int j = (i + 1) % V; //  The next vertex
vec3 color(1, 1, 1);
DrawLineSDL(screen, projectedVertices[i], projectedVertices[j], color);
}
}
*/

void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels) //Fixa till den här, inte original
{
	cout << "COMPUTE: " << vertexPixels.size() << " ms." << endl;
	// Find max and min y-value of the polygon
	// and compute the number of rows it occupies.
	int minimum = numeric_limits<int>::max();
	int maximum = numeric_limits<int>::min();
	for (int i = 0; i < vertexPixels.size(); i++) // Loop over them all
	{
		minimum = glm::min(vertexPixels[i].y, minimum); // Should be safer and faster
		maximum = glm::max(vertexPixels[i].y, maximum);
	}
	// Resize leftPixels and rightPixels
	// so that they have an element for each row.
	leftPixels.resize(maximum - minimum + 1);
	rightPixels.resize(maximum - minimum + 1);

	// Initialize the x-coordinates in leftPixels
	// to some really large value and the x-coordinates
	// in rightPixels to some really small value.
	for (int i = 0; i < leftPixels.size(); ++i)
	{
		leftPixels[i].x = SCREEN_WIDTH;
		leftPixels[i].y = i + minimum;

		rightPixels[i].x = 0;
		rightPixels[i].y = i  + minimum ;
	}

	// Loop through all edges of the polygon and use
	// linear interpolation to find the x-coordinate for
	// each row it occupies. Update the corresponding
	// values in rightPixels and leftPixels.
	//test
	for (int i = 0; i < vertexPixels.size(); ++i) { //skog

		int j = (i + 1) % vertexPixels.size(); // go to next vertex
		int pixels = abs(vertexPixels[i].y - vertexPixels[j].y) + 1;
		vector<Pixel> line(pixels);

		Interpolate(vertexPixels[i], vertexPixels[j], line);
		for (int k = 0; k < line.size(); ++k) {
			int rowNum = line[k].y - minimum;
			if (leftPixels[rowNum].x > line[k].x) {
				leftPixels[rowNum] = line[k];
			}
			if (rightPixels[rowNum].x < line[k].x) {
				rightPixels[rowNum] = line[k];

			}

		}
	}
}




void DrawRows(const vector<Pixel>& leftPixels, const vector<Pixel>& rightPixels) { //wibke
	for (int i = 0; i < leftPixels.size(); ++i) { //loop over rows
		DrawLineSDL(screen, leftPixels[i], rightPixels[i], currentColor);
	}
}

void DrawPolygon(const vector<Vertex>& vertices) //instruction
{
	int V = vertices.size();
	vector<Pixel> vertexPixels(V);
	for (int i = 0; i < V; ++i)
		VertexShader(vertices[i], vertexPixels[i]);

	vector<Pixel> leftPixels(3);
	vector<Pixel> rightPixels(3);
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	DrawRows(leftPixels, rightPixels);
}


int main(int argc, char* argv[])
{
	// Test case for ComputePolygonRows()
	/*
	vector<ivec2> vertexPixels(3);
	vertexPixels[0] = ivec2(10, 5);
	vertexPixels[1] = ivec2(5, 10);
	vertexPixels[2] = ivec2(15, 15);
	vector<ivec2> leftPixels;
	vector<ivec2> rightPixels;
	ComputePolygonRows(vertexPixels, leftPixels, rightPixels);
	for (size_t row = 0; row < leftPixels.size(); ++row) {
	cout << "Start: ("
	<< leftPixels[row].x << ","
	<< leftPixels[row].y << "). "
	<< "End: ("
	<< rightPixels[row].x << ","
	<< rightPixels[row].y << "). " << endl;
	}
	*/
	//start of main
	LoadTestModel(triangles);
	screen = InitializeSDL(SCREEN_WIDTH, SCREEN_HEIGHT);
	t = SDL_GetTicks();	// Set start value for timer.

	while (NoQuitMessageSDL())
	{
		Update();
		Draw();
	}

	SDL_SaveBMP(screen, "screenshot.bmp");
	return 0;
}

void Update()
{
	// Compute frame time:
	int t2 = SDL_GetTicks();
	float dt = float(t2 - t);
	t = t2;
	cout << "Render time: " << dt << " ms." << endl;

	Uint8* keystate = SDL_GetKeyState(0);

	R[0][0] = cos(yaw);
	R[0][1] = 0;
	R[0][2] = sin(yaw);

	//down
	R[1][0] = 0;
	R[1][1] = 1;
	R[1][2] = 0;

	//forward
	R[2][0] = -sin(yaw);
	R[2][1] = 0;
	R[2][2] = cos(yaw);

	float pan = 0.02f;

	if (keystate[SDLK_UP]) {
		cameraPos.z += pan;
	}
	if (keystate[SDLK_DOWN]) {
		cameraPos.z -= pan;
	}
	if (keystate[SDLK_LEFT]) {
		yaw += pan;
	}
	if (keystate[SDLK_RIGHT]) {
		yaw -= pan;
	}

	//Move lightsource
	if (keystate[SDLK_w])
		lightPos.y += pan;

	if (keystate[SDLK_s])
		lightPos.y -= pan;

	if (keystate[SDLK_d])
		lightPos.x += pan;

	if (keystate[SDLK_a])
		lightPos.x -= pan;

	if (keystate[SDLK_e])
		lightPos.z += pan;

	if (keystate[SDLK_q])
		lightPos.z -= pan;
}

void Draw()
{
	SDL_FillRect(screen, 0, 0);
	float apa = 1;
	if (SDL_MUSTLOCK(screen))
		SDL_LockSurface(screen);

	//clear depthbuffer
	for (int y = 0; y < SCREEN_HEIGHT; ++y) {
		for (int x = 0; x < SCREEN_WIDTH; ++x) {

			depthBuffer[y][x] = 0;
		}
	}

	//cout << "Start av drawaaaaaaaaaaaaaaaaaaa: " << triangles.size() << " ." << endl;
	for (int i = 0; i<triangles.size(); ++i)
	{
		currentColor = triangles[i].color;
		vector<Vertex> vertices(3); //vec3 -> Vertex step 6
		currentNormal = triangles[i].normal;
		currentReflectance = triangles[i].color;

		vertices[0].position = triangles[i].v0; //Adding .position to all the vertices in step 6								
		//cout << "Start av vertice0: " << triangles.size() << " ." << endl;
		vertices[1].position = triangles[i].v1;		
		//cout << "Start av vertice1: " << triangles.size() << " ." << endl;
		vertices[2].position = triangles[i].v2;		
		//cout << "Start av vertice2: " << triangles.size() << " ." << endl;
		//DrawPolygonEdges(vertices);
		DrawPolygon(vertices);
		//cout << "Slut av vertice2: " << triangles.size() << " ." << endl;
	}
	cout << "Start av vertice00000: " << triangles.size() << " ." << endl;
	if (SDL_MUSTLOCK(screen))
		SDL_UnlockSurface(screen);

	SDL_UpdateRect(screen, 0, 0, 0, 0);
	//cout << "Slut av draw: " << apa << " ." << endl;
}
