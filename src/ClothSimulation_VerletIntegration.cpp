#include <iostream>

#include <cmath>
#include <cstdint>
#include <SFML/Graphics/CircleShape.hpp>
#include <SFML/Graphics/Color.hpp>
#include <SFML/Graphics/RenderWindow.hpp>
#include <SFML/Graphics/Vertex.hpp>
#include <SFML/Window/Event.hpp>
#include <SFML/Window/VideoMode.hpp>

#include <array>
#include <SFML/Graphics.hpp>
#include <SFML/Graphics/PrimitiveType.hpp>
#include <SFML/Window/Mouse.hpp>

namespace
{
	const std::int32_t WIDTH = 2250;
	const std::int32_t HEIGHT = 1200;
	const float PARTICLE_RADIUS = 10.0f;
	const float GRAVITY = 15.0f;
	const float TIME_STEP = 0.1f;

	const int ROW = 10;
	const int COL = 10;
	const float REST_DISTANCE = 30.0f;

	const float CLICK_TOLERANCE = 5.0f;
}

class Particle
{
public:
	sf::Vector2f m_position;
	sf::Vector2f m_previous_position;
	sf::Vector2f m_acceleration;

	bool m_isPinned = false;

	Particle(float x, float y, bool pinned): m_position(x, y), m_previous_position(x, y), m_acceleration(0, 0), m_isPinned{pinned} {}

	auto applyForce(const sf::Vector2f &force) -> void
	{
		if (m_isPinned)
			return;

		m_acceleration += force;
	}

	auto update(float timeStep) -> void
	{
		if (m_isPinned)
			return;

		sf::Vector2f velocity_ = m_position - m_previous_position;
		m_previous_position = m_position;
		m_position += velocity_ + m_acceleration * timeStep * timeStep;
		m_acceleration = sf::Vector2f(0, 0);
	}

	auto constrainToBounds(float width, float height) -> void
	{
		if (m_position.x < 0)
			m_position.x = 0;

		if (m_position.x > width)
			m_position.x = width;

		if (m_position.y < 0)
			m_position.y = 0;

		if (m_position.y > height)
			m_position.y = height;
	}
};

class Constraint
{
public:
	Particle *m_p1;
	Particle *m_p2;
	float m_initialLength;
	bool m_active = true;

	Constraint(Particle *p1, Particle *p2): m_p1{p1}, m_p2{p2}, m_initialLength{}
	{
		m_initialLength = std::hypot(p2->m_position.x - p1->m_position.x, p2->m_position.y - p1->m_position.y);
	}

	auto satisfy() -> void
	{
		if (!m_active)
			return;

		sf::Vector2f delta_ = m_p2->m_position - m_p1->m_position;
		float currentLen_ = std::hypot(delta_.x, delta_.y);

		if (currentLen_ == 0.0f)
			currentLen_ = std::numeric_limits<float>::min();

		float diff_ = (currentLen_ - m_initialLength) / currentLen_;

		sf::Vector2f correction_ = delta_ * 0.5f * diff_;

		if (!m_p1->m_isPinned)
			m_p1->m_position += correction_;
		if (!m_p2->m_isPinned)
			m_p2->m_position -= correction_;
	}

	auto deactivate() -> void
	{
		m_active = false;
	}
};

class InputHandler
{
public:
	static auto handleMouseEvent(const std::optional<sf::Event> &event, std::vector<Particle> &particles, std::vector<Constraint> &constraints) -> void
	{
		if (const auto *mouseEvent_ = event->getIf<sf::Event::MouseButtonPressed>())
		{
			if (mouseEvent_->button == sf::Mouse::Button::Left)
			{
				auto mouse_x_ = static_cast<float>(mouseEvent_->position.x);
				auto mouse_y_ = static_cast<float>(mouseEvent_->position.y);
				tearCloth(mouse_x_, mouse_y_, particles, constraints);
			}
		}
	}

private:
	static auto point_to_segment_distance(float px, float py, float x1, float y1, float x2, float y2) -> float
	{
		float ABx = x2 - x1;
		float ABy = y2 - y1;

		float APx = px - x1;
		float APy = py - y1;

		float BPx = px - x2;
		float BPy = py - y2;

		float AB_AP = ABx * APx + ABy * APy;
		float AB_AB = ABx * ABx + ABy * ABy;
		float t = AB_AP / AB_AB;

		if (t < 0.0f)
		{
			return std::sqrt(APx * APx + APy * APy);
		}
		else if (t > 1.0f)
		{
			return std::sqrt(BPx * BPx + BPy * BPy);
		}
		else
		{
			float proj_x = x1 + t * ABx;
			float proj_y = y1 + t * ABy;
			return std::sqrt((px - proj_x) * (px - proj_x) + (py - proj_y) * (py - proj_y));
		}
	}

	static auto findNearestConstraint(float mouse_x, float mouse_y, const std::vector<Constraint> &constraints) -> Constraint *
	{
		Constraint *nearest_constraint_ = nullptr;
		float min_distance_ = CLICK_TOLERANCE;

		for (const auto &constraint : constraints)
		{
			float distance_ = point_to_segment_distance(mouse_x, mouse_y,
														constraint.m_p1->m_position.x, constraint.m_p1->m_position.y,
														constraint.m_p2->m_position.x, constraint.m_p2->m_position.y);
			if (distance_ < min_distance_)
			{
				min_distance_ = distance_;
				nearest_constraint_ = const_cast<Constraint *>(&constraint);
			}
		}
		return nearest_constraint_;
	}

	static auto tearCloth(float mouse_x, float mouse_y, const std::vector<Particle> &particles, std::vector<Constraint> &constraints) -> void
	{
		Constraint *nearest_ = findNearestConstraint(mouse_x, mouse_y, constraints);
		if (nearest_)
		{
			nearest_->deactivate();
		}
	}
};

int main()
{
	sf::RenderWindow window_{sf::VideoMode{sf::Vector2u{WIDTH, HEIGHT}}, "Cloth Simulation"};
	window_.setFramerateLimit(60);

	std::vector<Particle> particles_{};
	std::vector<Constraint> constraints_{};

	for (int row = 0; row < ROW; row++)
	{
		for (int col = 0; col < COL; col++)
		{
			float x = col * REST_DISTANCE + WIDTH / 3;
			float y = row * REST_DISTANCE + HEIGHT / 3;
			particles_.emplace_back(x, y, (row == 0));
		}
	}

	for (int row = 0; row < ROW; row++)
	{
		for (int col = 0; col < COL; col++)
		{
			if (col < COL - 1)
				constraints_.emplace_back(&particles_[row * COL + col], &particles_[row * COL + col + 1]);

			if (row < ROW - 1)
				constraints_.emplace_back(&particles_[row * COL + col], &particles_[(row + 1) * COL + col]);
		}
	}

	while (window_.isOpen())
	{
		while (const std::optional<sf::Event> event = window_.pollEvent())
		{
			if (event->is<sf::Event::Closed>())
				window_.close();

			InputHandler::handleMouseEvent(event, particles_, constraints_);
		}

		for (auto &particle_ : particles_)
		{
			particle_.applyForce(sf::Vector2f(0, GRAVITY));
			particle_.update(TIME_STEP);
			particle_.constrainToBounds(WIDTH, HEIGHT);
		}

		for (size_t i = 0; i < 1; i++)
		{
			for (auto &constraint_ : constraints_)
				constraint_.satisfy();
		}

		window_.clear(sf::Color::Black);

		//for (const auto &particle_ : particles_)
		//{
		//	sf::CircleShape circle_{PARTICLE_RADIUS};
		//	circle_.setFillColor(sf::Color::White);
		//	circle_.setPosition(particle_.m_position - sf::Vector2f(PARTICLE_RADIUS, PARTICLE_RADIUS));
		//	window_.draw(circle_);
		//}

		for (const auto &particle_ : particles_)
		{
			sf::Vertex point_{particle_.m_position, sf::Color::Red};
			window_.draw(&point_, 1, sf::PrimitiveType::Points);
		}

		for (const auto &constraint_ : constraints_)
		{
			if (!constraint_.m_active)
				continue;

			std::array<sf::Vertex, 2> line = {
				sf::Vertex{constraint_.m_p1->m_position, sf::Color::Red},
				sf::Vertex{constraint_.m_p2->m_position, sf::Color::Red}
			};

			window_.draw(line.data(), line.size(), sf::PrimitiveType::Lines);
		}

		window_.display();
	}
}
