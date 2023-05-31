/*

 Author: Yannis STEFANELLI

 Creation Date: 30-05-2023 10:48:03

 Description : Main backend rendering engine, uses SFML window & builtin OpenGL implementations

*/

#include <SFML/Graphics.hpp>
#include <iostream>
#include <string.h>

using namespace sf;

int main(int ac, char **av) {

    if (ac < 2) {
        std::cout << "Usage : ./shader-demo <shaderName>" << std::endl;
        return EXIT_SUCCESS;
    }

	const float winW = 1920;
	const float winH = 1080;
    
    ContextSettings settings;
    settings.antialiasingLevel = 8; // <---- Disable when you use heavy shaders

	RenderWindow window(VideoMode(winW, winH), "SFML Shader Example", Style::Close | Style::Fullscreen, settings);
	window.setMouseCursorVisible(false); // hide the cursor
    window.setFramerateLimit(60);
    window.setVerticalSyncEnabled(true);

	// Create a texture and a sprite for the shader
	Texture tex;
	tex.create(winW, winH);
	Sprite spr(tex);

	Shader shader;

    std::string shaderfile = "../shaders/";
    shaderfile += strcat(av[1], ".glsl");

	shader.loadFromFile(shaderfile, Shader::Fragment); // load the shader

	if (!shader.isAvailable()) {
		std::cout << "The shader is not available\n";
	}

	// Set the resolution parameter (the resoltion is divided to make the fire smaller)
	shader.setUniform("iResolution", Vector2f(winW / 2, winH / 2));

	// Use a timer to obtain the time elapsed
	Clock clk;
	clk.restart(); // start the timer
	
	sf::Vector2f fragCoord;

	while (window.isOpen()) {
		// Event handling
		Event event;

		while (window.pollEvent(event)) {
			// Exit the app when a key is pressed
			if (event.type == Event::KeyPressed) 
				return EXIT_SUCCESS;
		}

		// Set the others parameters who need to be updated every frames
		shader.setUniform("iTime", clk.getElapsedTime().asSeconds());

		Vector2i mousePos = Mouse::getPosition(window);
		shader.setUniform("iMouse", Vector2f(mousePos.x, mousePos.y - winH/2));

		// Draw the sprite with the shader on it
		window.clear();
		
		window.draw(spr, &shader);
		window.display();
	}

	return 0;
}